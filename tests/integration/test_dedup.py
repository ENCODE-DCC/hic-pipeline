from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_dedup")
def test_dedup_alignable_bam_match(workflow_dir, bam_md5_is_expected):
    bam_path = workflow_dir / Path("test-output/result_alignable_deup.bam")
    assert bam_md5_is_expected(bam_path, "a57b13a066a691628910d11f59e3cf0d")
