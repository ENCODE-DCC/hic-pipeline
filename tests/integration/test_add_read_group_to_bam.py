from pathlib import Path

import pysam
import pytest


@pytest.mark.workflow("test_add_read_group_to_bam")
def test_read_group_added_to_bam_header(workflow_dir):
    bam_path = workflow_dir / Path("test-output/with_read_group.bam")
    with pysam.AlignmentFile(str(bam_path)) as bam:
        assert bam.header.get("RG")[0]["ID"] == "foo"
