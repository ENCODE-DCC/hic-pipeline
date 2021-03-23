from pathlib import Path

import pytest


@pytest.mark.workflow("test_align")
def test_align_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/aligned.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "993e6ebf0dffa0f5e1bd9c974d4a6795"
