from pathlib import Path

import pytest


@pytest.mark.workflow("test_hic")
def test_hic_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "bedd4a49ac87457775089d68b070c865"


@pytest.mark.workflow("test_hic")
def test_hic_pairs_match(workflow_dir, skip_n_lines_md5):
    pairs_path = workflow_dir / Path("test-output/pairix.bsorted.pairs.gz")
    pairs_md5 = skip_n_lines_md5(pairs_path, n_lines=6)
    assert pairs_md5 == "0e82a8e6db855a237506306e299d0d5f"
