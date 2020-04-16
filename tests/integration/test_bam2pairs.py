from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_bam2pairs")
def test_bam2pairs_pairs_match(workflow_dir, skip_n_lines_md5):
    pairs_path = workflow_dir / Path("test-output/pairix.bsorted.pairs.gz")
    pairs_md5 = skip_n_lines_md5(pairs_path, n_lines=5)
    assert pairs_md5 == "454ac5c291c8b154e152530efa084f1b"
