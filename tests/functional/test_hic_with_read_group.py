from pathlib import Path

import pysam
import pytest


@pytest.mark.workflow("test_hic_with_read_group")
def test_read_group_added_to_bam_header(workflow_dir):
    bam_path = workflow_dir / Path("test-output/aligned.bam")
    with pysam.AlignmentFile(str(bam_path)) as bam:
        assert bam.header.get("RG")[0]["ID"] == "foo"


@pytest.mark.workflow("test_hic_with_read_group")
def test_hic_with_read_group_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/merged_dedup.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "4de9898216a3fc817850d2e1459d3ccf"
