from contextlib import suppress as does_not_raise

import pytest

from hic_pipeline.jsonify_stats import (
    clean_key,
    jsonify_stats,
    parse_to_dict,
    parse_to_int_or_float,
    parse_to_multiple_percentages,
    parse_to_total_and_percentage,
    parse_to_total_and_two_percentages,
)


def test_parse_to_dict():
    data = ["foo bar: baz\n", "Baz Qux: quuux\n", "corge:\n"]
    result = parse_to_dict(data)
    assert result == {"foo_bar": "baz", "baz_qux": "quuux"}


def test_parse_to_total_and_percentage():
    result = parse_to_total_and_percentage("26 (0.01%)")
    assert result == {"total": 26, "pct": 0.01}


def test_parse_to_total_and_two_percentages():
    result = parse_to_total_and_two_percentages(
        "6,969 (2.09% / 2.17%)",
        first_percentage_name="foo",
        second_percentage_name="bar",
    )
    assert result == {"total": 6969, "foo": 2.09, "bar": 2.17}


def test_parse_to_total_and_two_percentages_slash_separator():
    result = parse_to_total_and_two_percentages(
        "1,458 (92.96% / 96.24%)",
        first_percentage_name="foo",
        second_percentage_name="bar",
    )
    assert result == {"total": 1458, "foo": 92.96, "bar": 96.24}


def test_parse_to_multiple_percentages():
    result = parse_to_multiple_percentages(
        "25% - 23% - 27% - 25%", fields=["first", "second", "third", "fourth"]
    )
    assert result == {"first": 25.0, "second": 23.0, "third": 27.0, "fourth": 25.0}


@pytest.mark.parametrize(
    "value,condition,expected",
    [
        ("32", does_not_raise(), 32),
        ("32.0", does_not_raise(), 32.0),
        ("N/A", pytest.raises(ValueError), "unreachable"),
    ],
)
def test_parse_to_int_or_float(value, condition, expected):
    with condition:
        result = parse_to_int_or_float(value)
        assert result == expected


@pytest.mark.parametrize(
    "key,expected",
    [
        (" Ligation Motif Present", "ligation_motif_present"),
        ("3' Bias (Long Range)", "3_prime_bias_long_range"),
        ("<500BP", "short_range_less_than_500bp"),
        ("500BP-5kB", "short_range_500bp_to_5kb"),
        ("Long Range (>20Kb)", "long_range_greater_than_20kb"),
        ("Alignable (Normal+Chimeric Paired)", "alignable_normal_and_chimeric_paired"),
        ("Intra-fragment Reads", "intra_fragment_reads"),
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
        "single_alignment": "496 (0.15%)",
        "average_insert_size": "0.00",
        "alignable_normal_and_chimeric_paired": "321559 (96.60%)",
        "unique_reads": "321559 (100.00%, 96.60%)",
        "duplicates": "0 (0.00%, 0.00%)",
        "library_complexity_estimate": "N/A",
        "intra_fragment_reads": "6,969 (2.09% / 2.17%)",
        "below_mapq_threshold": "309,458 (92.96% / 96.24%)",
        "hic_contacts": "5,132 (1.54% / 1.60%)",
        # "ligation_motif_present": "96 (0.03%)",
        "ligation_motif_present": "6 (0.00% / 0.00%)",
        "3_prime_bias_long_range": "N/A",
        "pair_type_percent_lior": "25% - 23% - 27% - 25%",
        "lior_convergence": "10000000000",
        "inter_chromosomal": "6 (0.00% / 0.00%)",
        "intra_chromosomal": "5,126 (1.54% / 1.59%)",
        "short_range_less_than_500bp": "1,040 (0.31% / 0.32%)",
        "short_range_500bp_to_5kb": "3,387 (1.02% / 1.05%)",
        "short_range_5kb_to_20kb": "110 (0.03% / 0.03%)",
        "long_range_greater_than_20kb": "589 (0.18% / 0.18%)",
    }

    result = jsonify_stats(parsed_data)
    assert result == {
        "short_range_500bp_to_5kb": {
            "pct_of_sequenced_reads": 1.05,
            "pct_of_unique_reads": 1.02,
            "total": 3387,
        },
        "short_range_5kb_to_20kb": {
            "pct_of_sequenced_reads": 0.03,
            "pct_of_unique_reads": 0.03,
            "total": 110,
        },
        "alignable_(normal+chimeric_paired)": {"pct": 96.6, "total": 321559},
        "average_insert_size": 0.0,
        "below_mapq_threshold": {
            "pct_of_sequenced_reads": 96.24,
            "pct_of_unique_reads": 92.96,
            "total": 309458,
        },
        "chimeric_ambiguous": {"pct": 0.01, "total": 26},
        "chimeric_paired": {"pct": 0.0, "total": 1},
        "duplicates": {
            "pct_of_alignable_reads": 0.0,
            "pct_of_sequenced_reads": 0.0,
            "total": 0,
        },
        "hic_contacts": {
            "pct_of_sequenced_reads": 1.6,
            "pct_of_unique_reads": 1.54,
            "total": 5132,
        },
        "inter_chromosomal": {
            "pct_of_sequenced_reads": 0.0,
            "pct_of_unique_reads": 0.0,
            "total": 6,
        },
        "intra_chromosomal": {
            "pct_of_sequenced_reads": 1.59,
            "pct_of_unique_reads": 1.54,
            "total": 5126,
        },
        "intra_fragment_reads": {
            "pct_of_sequenced_reads": 2.17,
            "pct_of_unique_reads": 2.09,
            "total": 6969,
        },
        "short_range_less_than_500bp": {
            "pct_of_sequenced_reads": 0.32,
            "pct_of_unique_reads": 0.31,
            "total": 1040,
        },
        "ligation_motif_present": {
            "pct_of_sequenced_reads": 0.0,
            "pct_of_unique_reads": 0.0,
            "total": 6,
        },
        "lior_convergence": 10000000000,
        "long_range_greater_than_20kb": {
            "pct_of_sequenced_reads": 0.18,
            "pct_of_unique_reads": 0.18,
            "total": 589,
        },
        "normal_paired": {"pct": 96.6, "total": 321558},
        "pair_type_percent_lior": {
            "pct_inner": 23.0,
            "pct_left": 25.0,
            "pct_outer": 27.0,
            "pct_right": 25.0,
        },
        "sequenced_read_pairs": 332888,
        "single_alignment": {"pct": 0.15, "total": 496},
        "unique_reads": {
            "pct_of_alignable_reads": 100.0,
            "pct_of_sequenced_reads": 96.6,
            "total": 321559,
        },
        "unmapped": {"pct": 3.4, "total": 11303},
    }
