import pytest

from hic_pipeline.normalize_assembly_name import (
    get_wdl_boolean_string,
    normalize_assembly_name,
    normalize_grch_name,
)


@pytest.mark.parametrize(
    "assembly_name,expected",
    [
        ("GRCh38", ("hg38", True)),
        ("hg19", ("hg19", True)),
        ("GRCh37", ("hg19", True)),
        ("grch36", ("hg18", True)),
        ("unknown", ("unknown", False)),
    ],
)
def test_normalize_assembly_name(assembly_name, expected):
    result = normalize_assembly_name(assembly_name)
    assert result == expected


@pytest.mark.parametrize(
    "assembly_name,expected",
    [("GRCh38", "hg38"), ("GRCh37", "hg19"), ("GRCh36", "hg18")],
)
def test_normalize_grch_name_name(assembly_name, expected):
    result = normalize_grch_name(assembly_name)
    assert result == expected


@pytest.mark.parametrize("boolean,expected", [(True, "true"), (False, "false")])
def test_get_wdl_boolean_string(boolean, expected):
    result = get_wdl_boolean_string(boolean)
    assert result == expected
