from pathlib import Path

import pytest


@pytest.mark.workflow("test_ultima")
def test_ultima_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "59df19e09fc98126c3265984c07b0db5"
