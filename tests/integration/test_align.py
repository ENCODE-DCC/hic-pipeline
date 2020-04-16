from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_align")
def test_align_alignable_bam_match(workflow_dir, bam_md5_is_expected):
    bam_path = workflow_dir / Path("test-output/alignable.bam")
    assert bam_md5_is_expected(bam_path, "5682b857978a0e988523f0c979033d73")


@pytest.mark.workflow(name="test_align")
def test_align_mapq0_bam_match(workflow_dir, bam_md5_is_expected):
    bam_path = workflow_dir / Path("test-output/unmapped.bam")
    assert bam_md5_is_expected(bam_path, "e24ff1e0d5e6471fec6976d0aeaa1061")
