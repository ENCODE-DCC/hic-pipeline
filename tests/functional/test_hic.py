from pathlib import Path

import pytest


@pytest.mark.workflow("test_hic")
def test_hic_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/result_alignable_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "0bc6522586908f128d0f24a4dafe4065"


@pytest.mark.workflow("test_hic")
def test_hic_pairs_match(workflow_dir, skip_n_lines_md5):
    pairs_path = workflow_dir / Path("test-output/pairix.bsorted.pairs.gz")
    pairs_md5 = skip_n_lines_md5(pairs_path, n_lines=5)
    assert pairs_md5 == "a95ca169256c73b02f2c347b85b41215"
