import textwrap
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
    process_data,
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
        ("Average Insert Size", "avg_insert_size"),
        ("Unmapped", "unmapped_reads"),
        ("One or both reads unmapped", "one_or_both_reads_unmapped"),
        ("2 alignments", "2_alignments"),
        (" 2 alignments (A...B)", "2_alignments_a_b"),
        (" 2 alignments (A1….A2B; A1B2...B1A2)", "2_alignments_a1_a2b_a1b2_b1a2"),
        ("3 or more alignments", "3_or_more_alignments"),
        (
            "Library Complexity Estimate (1 alignment)*",
            "library_complexity_estimate_1_alignment",
        ),
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
        "unmapped_reads": "11303 (3.40%)",
        "single_alignment": "496 (0.15%)",
        "avg_insert_size": "0.00",
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
        "avg_insert_size": 0.0,
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
        "pct_unmapped_reads": 3.4,
        "unmapped_reads": 11303,
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


def test_process_paired_ended():
    data = textwrap.dedent(
        """Read type: Paired End
        Sequenced Read Pairs:  4079509
        No chimera found: 8979 (0.22%)
         One or both reads unmapped: 8979 (0.22%)
        2 alignments: 2522618 (61.84%)
         2 alignments (A...B): 0 (0.00%)
         2 alignments (A1….A2B; A1B2...B1A2): 2522618 (61.84%)
        3 or more alignments: 31489 (0.77%)
        Ligation Motif Present: N/A
        Average insert size: 0.00
        Total Unique: 2473835 (98.07%, 60.64%)
        Total Duplicates: 48783 (1.93%, 1.20%)
        Library Complexity Estimate*: 64,379,945
        Intra-fragment Reads: N/A
        Below MAPQ Threshold: 219,330 (5.38% / 8.87%)
        Hi-C Contacts: 2,254,505 (55.26% / 91.13%)
         3' Bias (Long Range): N/A
         Pair Type %(L-I-O-R): 25% - 25% - 25% - 25%
         L-I-O-R Convergence: 8111
        Inter-chromosomal: 511,681 (12.54% / 20.68%)
        Intra-chromosomal: 1,742,824 (42.72% / 70.45%)
        Short Range (<20Kb):
          <500BP: 220,276 (5.40% / 8.90%)
          500BP-5kB: 484,303 (11.87% / 19.58%)
          5kB-20kB: 173,245 (4.25% / 7.00%)
        Long Range (>20Kb): 865,000 (21.20% / 34.97%)
        """
    )
    result = process_data(data)
    assert result == {
        "run_type": "paired-ended",
        "sequenced_read_pairs": 4079509,
        "no_chimera_found": 8979,
        "pct_no_chimera_found": 0.22,
        "one_or_both_reads_unmapped": 8979,
        "pct_one_or_both_reads_unmapped": 0.22,
        "2_alignments": 2522618,
        "pct_2_alignments": 61.84,
        "2_alignments_a_b": 0,
        "pct_2_alignments_a_b": 0.0,
        "2_alignments_a1_a2b_a1b2_b1a2": 2522618,
        "pct_2_alignments_a1_a2b_a1b2_b1a2": 61.84,
        "3_or_more_alignments": 31489,
        "pct_3_or_more_alignments": 0.77,
        "avg_insert_size": 0.0,
        "total_unique": 2473835,
        "pct_unique_total_unique": 98.07,
        "pct_sequenced_total_unique": 60.64,
        "total_duplicates": 48783,
        "pct_unique_total_duplicates": 1.93,
        "pct_sequenced_total_duplicates": 1.2,
        "library_complexity_estimate": 64379945,
        "below_mapq_threshold": 219330,
        "pct_unique_below_mapq_threshold": 5.38,
        "pct_sequenced_below_mapq_threshold": 8.87,
        "hic_contacts": 2254505,
        "pct_unique_hic_contacts": 55.26,
        "pct_sequenced_hic_contacts": 91.13,
        "pct_left_pair_type": 25.0,
        "pct_inner_pair_type": 25.0,
        "pct_outer_pair_type": 25.0,
        "pct_right_pair_type": 25.0,
        "lior_convergence": 8111,
        "inter_chromosomal": 511681,
        "pct_unique_inter_chromosomal": 12.54,
        "pct_sequenced_inter_chromosomal": 20.68,
        "intra_chromosomal": 1742824,
        "pct_unique_intra_chromosomal": 42.72,
        "pct_sequenced_intra_chromosomal": 70.45,
        "short_range_less_than_500bp": 220276,
        "pct_unique_short_range_less_than_500bp": 5.4,
        "pct_sequenced_short_range_less_than_500bp": 8.9,
        "short_range_500bp_to_5kb": 484303,
        "pct_unique_short_range_500bp_to_5kb": 11.87,
        "pct_sequenced_short_range_500bp_to_5kb": 19.58,
        "short_range_5kb_to_20kb": 173245,
        "pct_unique_short_range_5kb_to_20kb": 4.25,
        "pct_sequenced_short_range_5kb_to_20kb": 7.0,
        "long_range_greater_than_20kb": 865000,
        "pct_unique_long_range_greater_than_20kb": 21.2,
        "pct_sequenced_long_range_greater_than_20kb": 34.97,
    }


def test_process_single_ended():
    data = textwrap.dedent(
        """Read type: Single End
        Sequenced Reads:  4079509
        No chimera found: 1525402 (37.39%)
         0 alignments: 8979 (0.22%)
         1 alignment: 1516423 (37.17%)
        2 alignments: 2522618 (61.84%)
        3 or more alignments: 31489 (0.77%)
        Ligation Motif Present: N/A
        1 alignment unique: 1496596 (98.69%, 36.69%)
        1 alignment duplicates: 19827 (1.31%, 0.49%)
        2 alignment unique: 2493662 (98.85%, 61.13%)
        2 alignment duplicates: 28956 (1.15%, 0.71%)
        Total Unique: 3990258 (98.79%, 97.81%)
        Total Duplicates: 48783 (1.21%, 1.20%)
        Library Complexity Estimate (1 alignment)*: 57,483,498
        Library Complexity Estimate (2 alignments)*: 109,041,497
        Library Complexity Estimate (1+2 above)*: 166,524,995
        Intra-fragment Reads: N/A
        Below MAPQ Threshold: 381,946 (9.36% / 15.32%)
        Hi-C Contacts: 2,111,716 (51.76% / 84.68%)
         3' Bias (Long Range): N/A
         Pair Type %(L-I-O-R): 25% - 25% - 25% - 25%
         L-I-O-R Convergence: 8111
        Inter-chromosomal: 463,914 (11.37% / 18.60%)
        Intra-chromosomal: 1,647,802 (40.39% / 66.08%)
        Short Range (<20Kb):
          <500BP: 209,033 (5.12% / 8.38%)
          500BP-5kB: 459,907 (11.27% / 18.44%)
          5kB-20kB: 163,690 (4.01% / 6.56%)
        Long Range (>20Kb): 815,172 (19.98% / 32.69%)
        """
    )
    result = process_data(data)
    assert result == {
        "run_type": "single-ended",
        "sequenced_reads": 4079509,
        "no_chimera_found": 1525402,
        "pct_no_chimera_found": 37.39,
        "0_alignments": 8979,
        "pct_0_alignments": 0.22,
        "1_alignment": 1516423,
        "pct_1_alignment": 37.17,
        "2_alignments": 2522618,
        "pct_2_alignments": 61.84,
        "3_or_more_alignments": 31489,
        "pct_3_or_more_alignments": 0.77,
        "1_alignment_unique": 1496596,
        "pct_1_alignment_unique": 98.69,
        "pct_sequenced_1_alignment_unique": 36.69,
        "1_alignment_duplicates": 19827,
        "pct_1_alignment_duplicates": 1.31,
        "pct_sequenced_1_alignment_duplicates": 0.49,
        "2_alignment_unique": 2493662,
        "pct_2_alignment_unique": 98.85,
        "pct_sequenced_2_alignment_unique": 61.13,
        "2_alignment_duplicates": 28956,
        "pct_2_alignment_duplicates": 1.15,
        "pct_sequenced_2_alignment_duplicates": 0.71,
        "total_unique": 3990258,
        "pct_unique_total_unique": 98.79,
        "pct_sequenced_total_unique": 97.81,
        "total_duplicates": 48783,
        "pct_unique_total_duplicates": 1.21,
        "pct_sequenced_total_duplicates": 1.2,
        "library_complexity_estimate_1_alignment": 57483498,
        "library_complexity_estimate_2_alignments": 109041497,
        "library_complexity_estimate_1_and_2_alignments": 166524995,
        "below_mapq_threshold": 381946,
        "pct_unique_below_mapq_threshold": 9.36,
        "pct_sequenced_below_mapq_threshold": 15.32,
        "hic_contacts": 2111716,
        "pct_unique_hic_contacts": 51.76,
        "pct_sequenced_hic_contacts": 84.68,
        "pct_left_pair_type": 25.0,
        "pct_inner_pair_type": 25.0,
        "pct_outer_pair_type": 25.0,
        "pct_right_pair_type": 25.0,
        "lior_convergence": 8111,
        "inter_chromosomal": 463914,
        "pct_unique_inter_chromosomal": 11.37,
        "pct_sequenced_inter_chromosomal": 18.6,
        "intra_chromosomal": 1647802,
        "pct_unique_intra_chromosomal": 40.39,
        "pct_sequenced_intra_chromosomal": 66.08,
        "short_range_less_than_500bp": 209033,
        "pct_unique_short_range_less_than_500bp": 5.12,
        "pct_sequenced_short_range_less_than_500bp": 8.38,
        "short_range_500bp_to_5kb": 459907,
        "pct_unique_short_range_500bp_to_5kb": 11.27,
        "pct_sequenced_short_range_500bp_to_5kb": 18.44,
        "short_range_5kb_to_20kb": 163690,
        "pct_unique_short_range_5kb_to_20kb": 4.01,
        "pct_sequenced_short_range_5kb_to_20kb": 6.56,
        "long_range_greater_than_20kb": 815172,
        "pct_unique_long_range_greater_than_20kb": 19.98,
        "pct_sequenced_long_range_greater_than_20kb": 32.69,
    }
