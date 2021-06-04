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
    update_with_total_and_percentages,
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


def test_update_with_total_and_percentages():
    parsed = {"total": 1, "pct_foo": 3.0, "pct_bar": 5.0}
    output = {"baz": 3}
    parent_key = "qux"
    update_with_total_and_percentages(output, parsed, parent_key)
    assert output == {"baz": 3, "qux": 1, "pct_foo_qux": 3.0, "pct_bar_qux": 5.0}


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
        "pct_sequenced_short_range_500bp_to_5kb": 1.05,
        "pct_unique_short_range_500bp_to_5kb": 1.02,
        "short_range_500bp_to_5kb": 3387,
        "pct_sequenced_short_range_5kb_to_20kb": 0.03,
        "pct_unique_short_range_5kb_to_20kb": 0.03,
        "short_range_5kb_to_20kb": 110,
        "pct_alignable_normal_and_chimeric_paired": 96.6,
        "alignable_normal_and_chimeric_paired": 321559,
        "average_insert_size": 0.0,
        "pct_sequenced_below_mapq_threshold": 96.24,
        "pct_unique_below_mapq_threshold": 92.96,
        "below_mapq_threshold": 309458,
        "pct_chimeric_ambiguous": 0.01,
        "chimeric_ambiguous": 26,
        "pct_chimeric_paired": 0.0,
        "chimeric_paired": 1,
        "pct_alignable_duplicates": 0.0,
        "pct_sequenced_duplicates": 0.0,
        "duplicates": 0,
        "pct_sequenced_hic_contacts": 1.6,
        "pct_unique_hic_contacts": 1.54,
        "hic_contacts": 5132,
        "pct_sequenced_inter_chromosomal": 0.0,
        "pct_unique_inter_chromosomal": 0.0,
        "inter_chromosomal": 6,
        "pct_sequenced_intra_chromosomal": 1.59,
        "pct_unique_intra_chromosomal": 1.54,
        "intra_chromosomal": 5126,
        "pct_sequenced_intra_fragment_reads": 2.17,
        "pct_unique_intra_fragment_reads": 2.09,
        "intra_fragment_reads": 6969,
        "pct_sequenced_short_range_less_than_500bp": 0.32,
        "pct_unique_short_range_less_than_500bp": 0.31,
        "short_range_less_than_500bp": 1040,
        "pct_sequenced_ligation_motif_present": 0.0,
        "pct_unique_ligation_motif_present": 0.0,
        "ligation_motif_present": 6,
        "lior_convergence": 10000000000,
        "pct_sequenced_long_range_greater_than_20kb": 0.18,
        "pct_unique_long_range_greater_than_20kb": 0.18,
        "long_range_greater_than_20kb": 589,
        "pct_normal_paired": 96.6,
        "normal_paired": 321558,
        "pct_inner_pair_type": 23.0,
        "pct_left_pair_type": 25.0,
        "pct_outer_pair_type": 27.0,
        "pct_right_pair_type": 25.0,
        "sequenced_read_pairs": 332888,
        "pct_single_alignment": 0.15,
        "single_alignment": 496,
        "pct_alignable_unique_reads": 100.0,
        "pct_sequenced_unique_reads": 96.6,
        "unique_reads": 321559,
        "pct_unmapped": 3.4,
        "unmapped": 11303,
    }


def test_jsonify_stats_without_nas():
    parsed_data = {
        "ligation_motif_present": "1077293854 (49.22%)",
        "library_complexity_estimate": "9,171,866,565",
        "3_prime_bias_long_range": "96% - 4%",
    }

    result = jsonify_stats(parsed_data)
    assert result == {
        "pct_3_prime_bias_long_range": 96.0,
        "pct_5_prime_bias_long_range": 4.0,
        "library_complexity_estimate": 9171866565,
        "ligation_motif_present": 1077293854,
        "pct_ligation_motif_present": 49.22,
    }
