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

ALLOWED_STATUSES = ("released", "in progress")


def main():
    parser = get_parser()
    args = parser.parse_args()
    auth = read_auth_from_file(args.keypair_file)
    experiment = get_experiment(args.accession, auth=auth)
    fastqs = get_fastqs_from_experiment(experiment)
    input_json = get_input_json(
        fastqs=fastqs, assembly_name=args.assembly_name, enzyme=args.enzyme
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


def get_input_json(fastqs, assembly_name, enzyme):
    input_json = {
        "hic.fastq": fastqs,
        "hic.assembly_name": assembly_name,
        "hic.chrsz": REFERENCE_FILES[assembly_name]["chrom_sizes"],
        "hic.reference_index": REFERENCE_FILES[assembly_name]["bwa_index"],
        "hic.restriction_enzymes": [enzyme],
    }
    if enzyme != "none":
        input_json["hic.restriction_sites"] = REFERENCE_FILES[assembly_name][
            "restriction_sites"
        ][enzyme]
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
    parser.add_argument(
        "-e",
        "--enzyme",
        choices=("HindIII", "DpnII", "MboI", "none"),
        required=True,
        help="Name of restriction enzyme",
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
    return parser


if __name__ == "__main__":
    main()
