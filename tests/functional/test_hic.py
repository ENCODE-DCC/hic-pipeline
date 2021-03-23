from pathlib import Path

import pytest


@pytest.mark.workflow("test_hic")
def test_hic_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "6df22b9574332f4db46a03ee15fbd73a"


@pytest.mark.workflow("test_hic")
def test_hic_pairs_match(workflow_dir, skip_n_lines_md5):
    pairs_path = workflow_dir / Path("test-output/pairix.bsorted.pairs.gz")
    pairs_md5 = skip_n_lines_md5(pairs_path, n_lines=5)
    assert pairs_md5 == "16188e32caf5ee02fb98a6cb417af16d"
