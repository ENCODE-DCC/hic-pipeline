import hashlib

import pytest


@pytest.mark.workflow("test_multiple_libraries")
def test_multiple_libraries_hic_match(workflow_dir):
    hic_path = next(
        workflow_dir.glob(
            "hic/*/call-create_hic_with_chrom_sizes/shard-1/execution/inter_30.hic"
        )
    )
    hic_md5sum = hashlib.md5(hic_path.read_bytes()).hexdigest()
    assert hic_md5sum == "69901660d2245767bc8048c0738d26bb"
