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
    }
}
ENZYMES = ("HindIII", "DpnII", "MboI", "none")
ALLOWED_STATUSES = ("released", "in progress")


def main():
    parser = get_parser()
    args = parser.parse_args()
    auth = read_auth_from_file(args.keypair_file)
    experiment = get_experiment(args.accession, auth=auth)
    fastqs = get_fastqs_from_experiment(experiment)
    if args.ligation_site_regex is None:
        enzymes = get_enzymes_from_experiment(experiment)
        input_json = get_input_json(
            fastqs=fastqs, assembly_name=args.assembly_name, enzymes=enzymes
        )
    else:
        input_json = get_input_json(
            fastqs=fastqs,
            assembly_name=args.assembly_name,
            ligation_site_regex=args.ligation_site_regex,
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
    if len(fragmentation_methods) > 1:
        raise ValueError(
            "Currently only experiments with one fragmentation method are supported"
        )
    for enzyme in enzymes:
        if enzyme in fragmentation_methods[0]:
            used_enzymes.append(enzyme)
            break
    if not used_enzymes:
        raise ValueError(
            "Unsupported fragmentation method: {}".format(fragmentation_methods[0])
        )
    return used_enzymes


def get_fastqs_from_experiment(experiment):
    fastq_pairs_by_replicate = {}
    for file in experiment["files"]:
        if file["file_format"] == "fastq" and file["status"] in ALLOWED_STATUSES:
            biological_replicate = file["biological_replicates"][0]
            paired_with_id = file["paired_with"]
            paired_with_file = [
                f for f in experiment["files"] if f["@id"] == paired_with_id
            ][0]
            replicate_fastqs = fastq_pairs_by_replicate.get(biological_replicate)
            if replicate_fastqs is None:
                fastq_pairs_by_replicate[biological_replicate] = []
            if file["paired_end"] == "2":
                file, paired_with_file = paired_with_file, file
            fastq_pair = {
                "read_1": urljoin(PORTAL_URL, file["href"]),
                "read_2": urljoin(PORTAL_URL, paired_with_file["href"]),
            }
            if fastq_pair not in fastq_pairs_by_replicate[biological_replicate]:
                fastq_pairs_by_replicate[biological_replicate].append(fastq_pair)
    output = [replicate for replicate in fastq_pairs_by_replicate.values()]
    return output


def get_input_json(fastqs, assembly_name, enzymes=None, ligation_site_regex=None):
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
    parser.add_argument(
        "--keypair-file", help="Path to keypairs.json", default="~/keypairs.json"
    )
    parser.add_argument(
        "--assembly-name",
        choices=("GRCh38",),
        default="GRCh38",
        help="Name of assembly, mm10 is not yet supported",
    )
    parser.add_argument("--ligation-site-regex", help="Regex for ligation site")
    return parser


if __name__ == "__main__":
    main()
