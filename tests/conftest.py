import pytest

from pathlib import Path


@pytest.fixture
def test_data_dir():
    return Path('tests/data')


@pytest.fixture
def bams_match():
    from comparisons import compare_bams_as_sams
    return compare_bams_as_sams


@pytest.fixture
def skip_lines_match():
    from comparisons import skip_n_lines_and_compare
    return skip_n_lines_and_compare
