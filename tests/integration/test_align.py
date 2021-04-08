from pathlib import Path

import pysam
import pytest


@pytest.fixture
def bam_path(workflow_dir):
    return workflow_dir / Path("test-output/aligned.bam")


@pytest.mark.workflow("test_align")
def test_align_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/aligned.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "993e6ebf0dffa0f5e1bd9c974d4a6795"


@pytest.mark.workflow("test_align_with_read_group")
def test_align_with_read_group_read_group_added_to_bam_header(bam_path):
    with pysam.AlignmentFile(str(bam_path)) as bam:
        assert bam.header.get("RG")[0]["ID"] == "foo"
        assert bam.header.get("RG")[0]["DS"] == "cool"


@pytest.mark.workflow("test_align_with_read_group")
def test_align_with_read_group_bam_match(bam_path, bam_md5):
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "8cb6f293859775fdee81a51c73f51c1e"
