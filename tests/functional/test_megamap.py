from pathlib import Path

import pytest


@pytest.mark.workflow("test_megamap")
def test_megamap_snp_vcfs_match(workflow_dir, skip_n_lines_md5):
    vcf_path = workflow_dir / Path("test-output/snp.out.vcf.gz")
    assert skip_n_lines_md5(vcf_path, n_lines=230) == "706f8cd7abe84b915aeec5aaaa07d067"


@pytest.mark.workflow("test_megamap")
def test_megamap_indel_vcfs_match(workflow_dir, skip_n_lines_md5):
    vcf_path = workflow_dir / Path("test-output/indel.out.vcf.gz")
    assert skip_n_lines_md5(vcf_path, n_lines=228) == "a6caf602a322a61a2c61b62ce7c449bd"
