from pathlib import Path

import pytest


@pytest.mark.workflow("test_merge")
def test_merge_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "518937c3a4bf0111549af2d3fe9b3289"
