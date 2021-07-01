from pathlib import Path

import pytest


@pytest.mark.workflow("test_pre_to_pairs")
def test_bam2pairs_pairs_match(workflow_dir, skip_n_lines_md5):
    pairs_path = workflow_dir / Path("test-output/pairix.bsorted.pairs.gz")
    pairs_md5 = skip_n_lines_md5(pairs_path, n_lines=6)
    assert pairs_md5 == "771a978c2737820c9d8dee03e6f89655"
