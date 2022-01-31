import textwrap
from io import StringIO

from hic_pipeline.filter_chrom_sizes import filter_chrom_sizes


def test_filter_chrom_sizes():
    chrom_sizes_data = textwrap.dedent(
        """\
            chr1\t248956422
            chrM\t16569
            chr1_GL383518v1_alt\t182439
            chr1_KI270707v1_random\t32032
            chrUn_GL000195v1\t182896
            chrEBV\t171823
        """
    )
    chrom_sizes = StringIO(initial_value=chrom_sizes_data)
    output_file = StringIO()
    filter_chrom_sizes(chrom_sizes, output_file)
    assert output_file.getvalue() == "chr1\t248956422\nchrM\t16569\n"
