import argparse
import json
from pathlib import Path
from urllib.parse import urljoin

import requests

_PORTAL_URL = "https://www.encodeproject.org"
_ALLOWED_STATUSES = ("released", "in progress")
_REFERENCE_FILES = {
    "reference_fasta": urljoin(
        _PORTAL_URL,
        "/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
    ),
    "gatk_bundle_tar": "https://www.encodeproject.org/files/ENCFF975DYE/@@download/ENCFF975DYE.tar.gz",
}


def main():
    parser = _get_parser()
    args = parser.parse_args()
    auth = _read_auth_from_file(args.keypair_file)
    bams = _get_bams_from_experiment(
        accessions=args.accessions,
        auth=auth,
        use_submitted_file_names=args.use_submitted_file_names,
    )
    donor_id = _get_donor_id_from_experiment(
        accessions=args.accessions,
        auth=auth
    )
    input_json = _get_input_json(bams=bams, donor_id=donor_id)
    _write_json_to_file(input_json, args.outfile)


def _get_experiment(accession, auth=None):
    response = requests.get(
        urljoin(_PORTAL_URL, accession),
        auth=auth,
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    return response.json()


def _get_bams_from_experiment(accessions, auth=None, use_submitted_file_names=False):
    bams = []
    for accession in accessions:
        experiment = _get_experiment(accession, auth=auth)
        analysis = None
        for a in experiment["analyses"]:
            if (
                a["status"] in _ALLOWED_STATUSES
                and a["lab"] == "/labs/encode-processing-pipeline/"
                and a["pipelines"][0] == "/pipelines/ENCPL839OAB/"
            ):
                analysis = a
                break
        if analysis is None:
            raise ValueError("Could not find uniform processing pipeline analysis")
        for file in experiment["files"]:
            if (
                file["status"] in _ALLOWED_STATUSES
                and file["@id"] in analysis["files"]
                and file["file_format"] == "bam"
            ):
                if use_submitted_file_names:
                    bams.append(file["submitted_file_name"])
                else:
                    bams.append(urljoin(_PORTAL_URL, file["href"]))
    return bams


def _get_donor_id_from_experiment(accessions, auth=None):
    donor_id = None
    for accession in accessions:
        experiment = _get_experiment(accession, auth=auth)
        if donor_id is None:
            donor_id = experiment["replicates"][0]["library"]["biosample"]["donor"]["accession"]
        else:
            if donor_id != experiment["replicates"][0]["library"]["biosample"]["donor"]["accession"]:
                raise ValueError("These datasets do not belong to the same donor.")
    return donor_id


def _get_input_json(bams, donor_id):
    input_json = {
        "genophase.bams": bams,
        "genophase.gatk_bundle_tar": _REFERENCE_FILES["gatk_bundle_tar"],
        "genophase.reference_fasta": _REFERENCE_FILES["reference_fasta"],
        "genophase.donor_id": donor_id
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
        help="Experiments to pool for genophasing",
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
