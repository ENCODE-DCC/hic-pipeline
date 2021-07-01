version 1.0

import "../../../hic.wdl" as hic

workflow test_pre_to_pairs {
    input {
        File pre
        File chrom_sizes
    }

    call hic.pre_to_pairs { input:
        pre = pre,
        chrom_sizes = chrom_sizes,
    }
}
