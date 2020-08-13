import pytest

from hic_pipeline.get_ligation_site_regex import get_ligation_site_regex, get_parser


@pytest.mark.parametrize(
    "enzymes,expected",
    [
        (["HindIII"], "AAGCTAGCTT"),
        (["HindIII", "DpnII"], "(AAGCTAGCTT|GATCGATC)"),
        (["MboI", "DpnII"], "GATCGATC"),
    ],
)
def test_get_ligation_site_regex(enzymes, expected):
    result = get_ligation_site_regex(enzymes)
    assert result == expected


def test_parser():
    parser = get_parser()
    result = parser.parse_args(["-o", "foo", "--enzymes", "bar", "baz"])
    assert result.outfile == "foo"
    assert result.enzymes == ["bar", "baz"]
