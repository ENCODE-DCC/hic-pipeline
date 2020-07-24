from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_align")
def test_align_alignable_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/alignable.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "5682b857978a0e988523f0c979033d73"


@pytest.mark.workflow(name="test_align")
def test_align_unmapped_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/unmapped.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "e24ff1e0d5e6471fec6976d0aeaa1061"


@pytest.mark.workflow(name="test_align")
def test_align_mapq0_bam_match(workflow_dir, bam_md5):
    bam_path = workflow_dir / Path("test-output/mapq0.bam")
    bam_md5sum = bam_md5(bam_path)
    assert bam_md5sum == "e24ff1e0d5e6471fec6976d0aeaa1061"
