version 1.0

import "../../../hic.wdl" as hic

workflow test_create_accessibility_track {
    input {
        File pre
        File chrom_sizes
        Int num_cpus
    }

    call hic.create_accessiblity_track { input:
        pre = pre,
        chrom_sizes = chrom_sizes,
        num_cpus = num_cpus,
    }
}
