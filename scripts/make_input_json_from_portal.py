import argparse
import json
from pathlib import Path
from urllib.parse import urljoin

import requests

PORTAL_URL = "https://www.encodeproject.org"

REFERENCE_FILES = {
    "GRCh38": {
        "restriction_sites": {
            "HindIII": urljoin(
                PORTAL_URL, "/files/ENCFF984SUZ/@@download/ENCFF984SUZ.txt.gz"
            ),
            "DpnII": urljoin(
                PORTAL_URL, "/files/ENCFF132WAM/@@download/ENCFF132WAM.txt.gz"
            ),
            "MboI": urljoin(
                PORTAL_URL, "/files/ENCFF132WAM/@@download/ENCFF132WAM.txt.gz"
            ),
        },
        "bwa_index": urljoin(
            PORTAL_URL, "/files/ENCFF643CGH/@@download/ENCFF643CGH.tar.gz"
        ),
        "chrom_sizes": urljoin(
            PORTAL_URL,
            "/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv",
        ),
    },
    "mm10": {
        "restriction_sites": {
            "DpnII": urljoin(
                PORTAL_URL, "/files/ENCFF930KBK/@@download/ENCFF930KBK.txt.gz"
            ),
            "MboI": urljoin(
                PORTAL_URL, "/files/ENCFF930KBK/@@download/ENCFF930KBK.txt.gz"
            ),
        },
        "bwa_index": urljoin(
            PORTAL_URL, "/files/ENCFF018NEO/@@download/ENCFF018NEO.tar.gz"
        ),
        "chrom_sizes": urljoin(
            PORTAL_URL,
            "/files/mm10_no_alt.chrom.sizes/@@download/mm10_no_alt.chrom.sizes.tsv",
        ),
    },
}
ENZYMES = ("HindIII", "DpnII", "MboI", "MseI", "none")
_NO_ENZYME_FRAGMENTATION_METHODS = (
    "chemical (micrococcal nuclease)",
    "chemical (DNaseI)",
)
ALLOWED_STATUSES = ("released", "in progress")


def main():
    parser = get_parser()
    args = parser.parse_args()
    auth = read_auth_from_file(args.keypair_file)
    experiment = get_experiment(args.accession, auth=auth)
    organism = experiment["replicates"][0]["library"]["biosample"]["organism"]["@id"]
    if organism == "/organisms/mouse/":
        assembly_name = "mm10"
    elif organism == "/organisms/human/":
        assembly_name = "GRCh38"
    else:
        raise ValueError(f"Organism {organism} not supported")
    try:
        fastqs = get_fastqs_from_experiment(experiment)
    except ValueError:
        fastqs = get_fastqs_from_experiment(experiment, read_group_num_path_parts=2)
    if args.ligation_site_regex is None:
        enzymes = args.enzymes or get_enzymes_from_experiment(experiment)
        input_json = get_input_json(
            fastqs=fastqs,
            assembly_name=assembly_name,
            enzymes=enzymes,
            no_delta=args.no_delta,
            no_slice=args.no_slice,
        )
    else:
        input_json = get_input_json(
            fastqs=fastqs,
            assembly_name=assembly_name,
            ligation_site_regex=args.ligation_site_regex,
            no_delta=args.no_delta,
            no_slice=args.no_slice,
        )
    outfile = args.outfile or "{}.json".format(args.accession)
    write_json_to_file(input_json, outfile)


def get_experiment(accession, auth=None):
    response = requests.get(
        urljoin(PORTAL_URL, accession),
        auth=auth,
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    return response.json()


def get_enzymes_from_experiment(experiment, enzymes=ENZYMES):
    used_enzymes = []
    fragmentation_methods = []
    for replicate in experiment["replicates"]:
        fragmentation_methods.extend(replicate["library"]["fragmentation_methods"])
    fragmentation_methods = list(set(fragmentation_methods))
    for fragmentation_method in fragmentation_methods:
        if fragmentation_method in _NO_ENZYME_FRAGMENTATION_METHODS:
            used_enzymes.append("none")
        else:
            used_enzyme = [
                enzyme for enzyme in enzymes if enzyme in fragmentation_method
            ]
            if used_enzyme:
                used_enzymes += used_enzyme
            else:
                raise ValueError(
                    "Unsupported fragmentation method: {}".format(fragmentation_method)
                )
    return used_enzymes


def get_fastqs_from_experiment(experiment, read_group_num_path_parts=1):
    fastq_pairs_by_replicate = {}
    read_group_ids = set()
    for file in experiment["files"]:
        if file["file_format"] == "fastq" and file["status"] in ALLOWED_STATUSES:
            biological_replicate = file["biological_replicates"][0]
            paired_with = file.get("paired_with")
            # Same as how juicer.sh does it, assumes "_R1" is used to indicated read 1
            # It is perhaps a bit low level, but it is hard to construct a unique ID
            # per read pair from other pieces of portal metadata. fastq signature is
            # close but not readily useful
            library_accession = file["replicate"]["library"].split("/")[2]
            sample_name = experiment["accession"]
            read_group_id = "_".join(
                file["submitted_file_name"].split("/")[-read_group_num_path_parts:]
            )
            read_group_id = (
                read_group_id.replace(".fastq.gz", "")
                .replace("_R1", "")
                .replace("_R2", "")
            )
            if paired_with is not None:
                platform = "ILLUMINA"
                paired_with_file = [
                    f for f in experiment["files"] if f["@id"] == paired_with
                ][0]
                if file["paired_end"] == "2":
                    file, paired_with_file = paired_with_file, file
                fastq_pair = {
                    "read_1": urljoin(PORTAL_URL, file["href"]),
                    "read_2": urljoin(PORTAL_URL, paired_with_file["href"]),
                }
            else:
                platform = "LS454"
                fastq_pair = {
                    "read_1": urljoin(PORTAL_URL, file["href"]),
                }
            fastq_pair[
                "read_group"
            ] = f"@RG\\tID:{read_group_id}\\tSM:{sample_name}\\tPL:{platform}\\tLB:{library_accession}"
            replicate_fastqs = fastq_pairs_by_replicate.get(biological_replicate)
            if replicate_fastqs is None:
                fastq_pairs_by_replicate[biological_replicate] = []
            if fastq_pair not in fastq_pairs_by_replicate[biological_replicate]:
                if read_group_id not in read_group_ids:
                    read_group_ids.add(read_group_id)
                else:
                    raise ValueError(
                        f"Read group ids must be unique, found duplicate id {read_group_id}"
                    )
                fastq_pairs_by_replicate[biological_replicate].append(fastq_pair)
    output = [replicate for replicate in fastq_pairs_by_replicate.values()]
    return output


def get_input_json(
    fastqs,
    assembly_name,
    enzymes=None,
    ligation_site_regex=None,
    no_slice=False,
    no_delta=False,
):
    input_json = {
        "hic.fastq": fastqs,
        "hic.assembly_name": assembly_name,
        "hic.chrsz": REFERENCE_FILES[assembly_name]["chrom_sizes"],
        "hic.reference_index": REFERENCE_FILES[assembly_name]["bwa_index"],
    }
    if enzymes is not None:
        input_json["hic.restriction_enzymes"] = enzymes
        if enzymes != ["none"]:
            input_json["hic.restriction_sites"] = REFERENCE_FILES[assembly_name][
                "restriction_sites"
            ][enzymes[0]]

    if ligation_site_regex is not None:
        input_json["hic.ligation_site_regex"] = ligation_site_regex

    if "read_2" not in fastqs[0][0]:
        input_json["hic.intact"] = True

    if no_slice:
        input_json["hic.no_slice"] = True
    if no_delta:
        input_json["hic.no_delta"] = True

    return input_json


def write_json_to_file(data, outfile):
    Path(outfile).write_text(json.dumps(data, indent=2, sort_keys=True))


def read_auth_from_file(keypair_file):
    keypair_path = Path(keypair_file).expanduser()
    if keypair_path.exists():
        data = json.loads(keypair_path.read_text())
        return (data["submit"]["key"], data["submit"]["secret"])
    else:
        return None


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--accession",
        required=True,
        help="Accession of portal experiment to generate input for",
    )
    parser.add_argument("--outfile")
    parser.add_argument("--no-slice", action="store_true")
    parser.add_argument("--no-delta", action="store_true")
    parser.add_argument(
        "-e", "--enzymes", nargs="+", help="Restriction enzymes used in experiment"
    )
    parser.add_argument(
        "--keypair-file", help="Path to keypairs.json", default="~/keypairs.json"
    )
    parser.add_argument("--ligation-site-regex", help="Regex for ligation site")
    return parser


if __name__ == "__main__":
    main()
