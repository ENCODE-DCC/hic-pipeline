import argparse
import json
from pathlib import Path
from urllib.parse import urljoin

import requests

_PORTAL_URL = "https://www.encodeproject.org"
_ALLOWED_STATUSES = ("released", "in progress")
_REFERENCE_FILES = {
    "chrom_sizes": urljoin(
        _PORTAL_URL,
        "/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv",
    ),
}


def main():
    parser = _get_parser()
    args = parser.parse_args()
    auth = _read_auth_from_file(args.keypair_file)
    bigwig_files, hic_files = _get_bigwig_files_and_hic_files_from_experiment(
        accessions=args.accessions,
        auth=auth,
        use_submitted_file_names=args.use_submitted_file_names,
    )
    input_json = _get_input_json(bigwig_files=bigwig_files, hic_files=hic_files)
    _write_json_to_file(input_json, args.outfile)


def _get_experiment(accession, auth=None):
    response = requests.get(
        urljoin(_PORTAL_URL, accession),
        auth=auth,
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    return response.json()


def _get_bigwig_files_and_hic_files_from_experiment(
    accessions, auth=None, use_submitted_file_names=False
):
    bigwig_files = []
    hic_files = []
    for accession in accessions:
        experiment = _get_experiment(accession, auth=auth)
        analysis = None
        for a in experiment["analyses"]:
            if (
                a["status"] in _ALLOWED_STATUSES
                and a["lab"] == "/labs/encode-processing-pipeline/"
            ):
                analysis = a
                break
        if analysis is None:
            raise ValueError("Could not find uniform processing pipeline analysis")
        for file in experiment["files"]:
            if file["status"] in _ALLOWED_STATUSES and file["@id"] in analysis["files"]:
                if file["output_type"] == "DNA accessibility raw signal":
                    if use_submitted_file_names:
                        bigwig_files.append(file["submitted_file_name"])
                    else:
                        bigwig_files.append(urljoin(_PORTAL_URL, file["href"]))
                if (
                    file["output_type"]
                    == "mapping quality thresholded chromatin interactions"
                ):
                    if use_submitted_file_names:
                        hic_files.append(file["submitted_file_name"])
                    else:
                        hic_files.append(urljoin(_PORTAL_URL, file["href"]))
    return bigwig_files, hic_files


def _get_input_json(bigwig_files, hic_files):
    input_json = {
        "megamap.bigwig_files": bigwig_files,
        "megamap.hic_files": hic_files,
        "megamap.chrom_sizes": _REFERENCE_FILES["chrom_sizes"],
    }
    return input_json


def _write_json_to_file(data, outfile):
    Path(outfile).write_text(json.dumps(data, indent=2, sort_keys=True))


def _read_auth_from_file(keypair_file):
    keypair_path = Path(keypair_file).expanduser()
    if keypair_path.exists():
        data = json.loads(keypair_path.read_text())
        return (data["submit"]["key"], data["submit"]["secret"])
    else:
        return None


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile")
    parser.add_argument(
        "-a",
        "--accessions",
        nargs="+",
        help="Experiments to pool for megamap",
        required=True,
    )
    parser.add_argument(
        "--keypair-file", help="Path to keypairs.json", default="~/keypairs.json"
    )
    parser.add_argument(
        "-u",
        "--use-submitted-file-names",
        help="Use submitted file names (gs:// uris) for bams to avoid massive file download",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    main()
