from pathlib import Path

import pytest


@pytest.mark.workflow("test_hic")
def test_hic_merged_dedup_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "bedd4a49ac87457775089d68b070c865"
