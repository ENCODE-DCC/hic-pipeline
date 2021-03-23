from pathlib import Path

import pytest


@pytest.mark.workflow("test_dedup")
def test_dedup_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "58eaeef21ade4eabb71dfa0d3889a3be"
