import gzip
import hashlib
from pathlib import Path

import pysam
import pytest


@pytest.fixture
def test_data_dir():
    return Path("tests/data")


@pytest.fixture
def bam_md5():
    """
    Bams often don't match due to nondeterministic insertion of paths into the headers.
    Therefore we need to view it as headerless SAM for comparisons.
    """

    def _bam_md5(bam_path: Path) -> str:
        sam = pysam.view(str(bam_path))
        return md5sum(sam)

    return _bam_md5


@pytest.fixture
def skip_n_lines_md5():
    """
    Text files can sometimes contain nondeterministic data in the headers. This fixture
    returns a function that will compare the md5sums of a file after n lines have been
    skipped. Will decompress gzipped files if need be.
    """

    def _skip_n_lines_md5(file_path: Path, n_lines: int) -> str:
        try:
            with gzip.open(str(file_path), "rt") as f:
                lines = "\n".join(f.readlines()[n_lines:])
        except OSError:
            with open(str(file_path)) as f:
                lines = "\n".join(f.readlines()[n_lines:])
        return md5sum(lines)

    return _skip_n_lines_md5


def md5sum(file: str) -> str:
    return hashlib.md5(file.encode()).hexdigest()
