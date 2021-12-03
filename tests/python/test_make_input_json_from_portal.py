import pytest

from scripts.make_input_json_from_portal import (
    get_enzymes_from_experiment,
    get_fastqs_from_experiment,
    get_input_json,
)


def test_get_input_json():
    result = get_input_json(
        fastqs=[[{"read_1": "foo", "read_2": "bar"}]],
        assembly_name="GRCh38",
        enzymes=["MboI"],
    )
    assert result == {
        "hic.assembly_name": "GRCh38",
        "hic.chrsz": "https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv",
        "hic.fastq": [[{"read_1": "foo", "read_2": "bar"}]],
        "hic.reference_index": "https://www.encodeproject.org/files/ENCFF643CGH/@@download/ENCFF643CGH.tar.gz",
        "hic.restriction_enzymes": ["MboI"],
        "hic.restriction_sites": "https://www.encodeproject.org/files/ENCFF132WAM/@@download/ENCFF132WAM.txt.gz",
    }


def test_get_input_json_none_enzyme_has_no_restriction_sites():
    result = get_input_json(
        fastqs=[[{"read_1": "foo"}]], assembly_name="GRCh38", enzymes=["none"]
    )
    assert result == {
        "hic.assembly_name": "GRCh38",
        "hic.chrsz": "https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv",
        "hic.delta_resolutions": [1000, 5000, 10000],
        "hic.fastq": [[{"read_1": "foo"}]],
        "hic.reference_index": "https://www.encodeproject.org/files/ENCFF643CGH/@@download/ENCFF643CGH.tar.gz",
        "hic.restriction_enzymes": ["none"],
    }


def test_get_enzymes_from_experiment():
    result = get_enzymes_from_experiment(
        {
            "replicates": [
                {
                    "library": {
                        "fragmentation_methods": ["chemical (myenzyme restriction"]
                    }
                }
            ]
        },
        enzymes=["myenzyme"],
    )
    assert result == ["myenzyme"]


def test_get_enzymes_from_experiment_mnase_returns_none_enzyme():
    result = get_enzymes_from_experiment(
        {
            "replicates": [
                {
                    "library": {
                        "fragmentation_methods": ["chemical (micrococcal nuclease)"]
                    }
                }
            ]
        }
    )
    assert result == ["none"]


def test_get_enzymes_from_experiment_multiple_fragmentation_methods_raises():
    experiment = {
        "replicates": [
            {"library": {"fragmentation_methods": ["chemical (MboI restriction)"]}},
            {"library": {"fragmentation_methods": ["chemical (MseI restriction)"]}},
        ]
    }
    with pytest.raises(ValueError):
        get_enzymes_from_experiment(experiment, enzymes=["MboI", "MseI"])


def test_get_enzymes_from_experiment_unknown_fragmentation_methods_raises():
    experiment = {
        "replicates": [
            {"library": {"fragmentation_methods": ["chemical (myenzyme restriction"]}}
        ]
    }
    with pytest.raises(ValueError):
        get_enzymes_from_experiment(experiment, enzymes=["MboI", "MseI"])


def test_get_fastqs_from_experiment():
    experiment = {
        "accession": "ENCSR123ABC",
        "files": [
            {
                "@id": "foo",
                "biological_replicates": ["1"],
                "paired_end": "1",
                "paired_with": "bar",
                "status": "released",
                "file_format": "fastq",
                "href": "download1",
                "submitted_file_name": "/foo/bar_R1.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "bar",
                "biological_replicates": ["1"],
                "paired_end": "2",
                "paired_with": "foo",
                "status": "released",
                "file_format": "fastq",
                "href": "download2",
                "submitted_file_name": "/foo/bar_R2.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "baz",
                "biological_replicates": ["2"],
                "paired_end": "1",
                "paired_with": "qux",
                "status": "in progress",
                "file_format": "fastq",
                "href": "download3",
                "submitted_file_name": "/baz/qux_R1.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "qux",
                "biological_replicates": ["2"],
                "paired_end": "2",
                "paired_with": "baz",
                "status": "in progress",
                "file_format": "fastq",
                "href": "download4",
                "submitted_file_name": "/baz/qux_R2.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "quux",
                "biological_replicates": ["2"],
                "paired_end": "1",
                "paired_with": "corge",
                "status": "released",
                "file_format": "fastq",
                "href": "download5",
                "submitted_file_name": "/quux/corge_R1.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "corge",
                "biological_replicates": ["2"],
                "paired_end": "2",
                "paired_with": "quux",
                "status": "released",
                "file_format": "fastq",
                "href": "download6",
                "submitted_file_name": "/quux/corge_R2.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "grault",
                "biological_replicates": ["2"],
                "paired_end": "1",
                "paired_with": "garply",
                "status": "replaced",
                "file_format": "fastq",
                "href": "download7",
                "submitted_file_name": "/grault/garply_R1.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
            {
                "@id": "garply",
                "biological_replicates": ["2"],
                "paired_end": "2",
                "paired_with": "grault",
                "status": "replaced",
                "file_format": "fastq",
                "href": "download8",
                "submitted_file_name": "/grault/garply_R2.fastq.gz",
                "replicate": {
                    "library": "/libraries/foo/",
                },
            },
        ],
    }
    result = get_fastqs_from_experiment(experiment)
    assert result == [
        [
            {
                "read_1": "https://www.encodeproject.org/download1",
                "read_2": "https://www.encodeproject.org/download2",
                "read_group": "@RG\\tID:bar\\tSM:ENCSR123ABC\\tPL:ILLUMINA\\tLB:foo",
            }
        ],
        [
            {
                "read_1": "https://www.encodeproject.org/download3",
                "read_2": "https://www.encodeproject.org/download4",
                "read_group": "@RG\\tID:qux\\tSM:ENCSR123ABC\\tPL:ILLUMINA\\tLB:foo",
            },
            {
                "read_1": "https://www.encodeproject.org/download5",
                "read_2": "https://www.encodeproject.org/download6",
                "read_group": "@RG\\tID:corge\\tSM:ENCSR123ABC\\tPL:ILLUMINA\\tLB:foo",
            },
        ],
    ]
