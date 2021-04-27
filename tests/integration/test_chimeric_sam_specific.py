from pathlib import Path

import pytest


@pytest.mark.workflow("test_chimeric_sam_specific")
def test_dedup_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/chimeric_sam_specific.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "7f524c704e9fc1077423fdc09a813e26"
