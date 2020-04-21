from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_dedup")
def test_dedup_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/result_alignable_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "a57b13a066a691628910d11f59e3cf0d"
