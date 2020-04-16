from pathlib import Path

import pytest


@pytest.mark.workflow(name="test_merge")
def test_merge_bam_match(workflow_dir, bam_md5_is_expected):
    bam_path = workflow_dir / Path("test-output/merged_bam_files.bam")
    assert bam_md5_is_expected(bam_path, "5e0726f291740a645f514e7a13ed7d21")
