version 1.0

import "../../../hic.wdl" as hic

workflow test_pre_to_pairs {
    input {
        File pre
        File chrom_sizes
        RuntimeEnvironment runtime_environment
    }

    call hic.pre_to_pairs { input:
        pre = pre,
        chrom_sizes = chrom_sizes,
        runtime_environment = runtime_environment,
    }
}
