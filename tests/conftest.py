import gzip
import hashlib
from pathlib import Path

import pysam
import pytest


@pytest.fixture
def test_data_dir():
    return Path("tests/data")


@pytest.fixture
def bam_md5_is_expected():
    """
    Bams often don't match due to nondeterministic insertion of paths into the headers.
    Therefore we need to view it as headerless SAM for comparisons.
    """

    def _bam_md5_is_expected(bam_path: Path, expected_md5: str) -> bool:
        sam = pysam.view(str(bam_path))
        return md5sum(sam) == expected_md5

    return _bam_md5_is_expected


@pytest.fixture
def skip_n_lines_md5_is_expected():
    """
    Text files can sometimes contain nondeterministic data in the headers. This fixture
    returns a function that will compare the md5sums of a file after n lines have been
    skipped. Will decompress gzipped files if need be.
    """

    def _skip_n_lines_md5_is_expected(
        file_path: Path, expected_md5: str, n_lines: int
    ) -> bool:
        try:
            with gzip.open(str(file_path), "rt") as f:
                lines = "\n".join(f.readlines()[n_lines:])
        except OSError:
            with open(str(file_path)) as f:
                lines = "\n".join(f.readlines()[n_lines:])
        return md5sum(lines) == expected_md5

    return _skip_n_lines_md5_is_expected


def md5sum(file: str) -> str:
    return hashlib.md5(file.encode()).hexdigest()
