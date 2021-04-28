from pathlib import Path

import pytest


@pytest.mark.workflow("test_chimeric_sam_nonspecific")
def test_dedup_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/chimeric_sam_nonspecific.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "56d10d96daa12ec92542bc5d61ab3d72"
