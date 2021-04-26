import hashlib

import pytest


@pytest.mark.workflow("test_multiple_libraries")
def test_multiple_libraries_hic_match(workflow_dir):
    hic_path = next(
        workflow_dir.glob("hic/*/call-create_hic/shard-1/execution/inter_30.hic")
    )
    hic_md5sum = hashlib.md5(hic_path.read_bytes()).hexdigest()
    assert hic_md5sum == "d7186dc77cf8ff79c8efd4ba3ace1afb"
