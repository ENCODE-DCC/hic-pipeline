import pytest

from hic_pipeline.jsonify_stats import clean_key, jsonify_stats, parse_to_dict


def test_parse_to_dict():
    data = ["foo bar: baz\n", "Baz Qux: quuux\n", "corge:\n"]
    result = parse_to_dict(data)
    assert result == {"foo_bar": "baz", "baz_qux": "quuux"}


@pytest.mark.parametrize(
    "key,expected",
    [
        (" Ligation Motif Present", "ligation_motif_present"),
        ("3' Bias (Long Range)", "3_prime_bias_long_range"),
        ("<500BP", "less_than_500bp"),
        ("500BP-5kB", "500bp_to_5kb"),
        ("Long Range (>20Kb)", "long_range_greater_than_20kb"),
        ("Alignable (Normal+Chimeric Paired)", "alignable_normal_and_chimeric_paired"),
        ('Intra-fragment Reads', "intra_fragment_reads"),
        ("Hi-C Contacts", "hic_contacts"),
        ("Pair Type %(L-I-O-R)", "pair_type_percent_lior"),
    ],
)
def test_clean_key(key, expected):
    result = clean_key(key)
    assert result == expected


def test_jsonify_stats():
    parsed_data = {
        "sequenced_read_pairs": "332888",
        "normal_paired": "321558 (96.60%)",
        "chimeric_paired": "1 (0.00%)",
        "chimeric_ambiguous": "26 (0.01%)",
        "unmapped": "11303 (3.40%)",
        # "ligation_motif_present": "96 (0.03%)",
        "single_alignment": "496 (0.15%)",
        "average_insert_size": "0.00",
        "alignable_(normal+chimeric_paired)": "321559 (96.60%)",
        "unique_reads": "321559 (100.00%, 96.60%)",
        "duplicates": "0 (0.00%, 0.00%)",
        "library_complexity_estimate": "N/A",
        "intra-fragment_reads": "6,969 (2.09% / 2.17%)",
        "below_mapq_threshold": "309,458 (92.96% / 96.24%)",
        "hi-c_contacts": "5,132 (1.54% / 1.60%)",
        "ligation_motif_present": "6 (0.00% / 0.00%)",
        "3'_bias_(long_range)": "N/A",
        "pair_type_%(l-i-o-r)": "25% - 23% - 27% - 25%",
        "l-i-o-r_convergence": "10000000000",
        "inter-chromosomal": "6 (0.00% / 0.00%)",
        "intra-chromosomal": "5,126 (1.54% / 1.59%)",
        "<500bp": "1,040 (0.31% / 0.32%)",
        "500bp-5kb": "3,387 (1.02% / 1.05%)",
        "5kb-20kb": "110 (0.03% / 0.03%)",
        "long_range_(>20kb)": "589 (0.18% / 0.18%)",
    }

    result = jsonify_stats(parsed_data)
    assert result == {"foo": "bar"}
