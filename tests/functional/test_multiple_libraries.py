from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_multiple_libraries")
def test_multiple_libraries_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/result_alignable_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "1d36f1e13987f5de4df7cd08bcf7982d"


@pytest.mark.workflow(name="test_multiple_libraries")
def test_multiple_libraries_pairs_match(workflow_dir, skip_n_lines_md5):
    pairs_path = workflow_dir / Path("test-output/pairix.bsorted.pairs.gz")
    pairs_md5 = skip_n_lines_md5(pairs_path, n_lines=5)
    assert pairs_md5 == "d5896f7dd95a90bf35870e2c5592cb37"
