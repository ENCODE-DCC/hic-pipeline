version 1.0

import "../../../hic.wdl" as hic

workflow test_create_accessibility_track {
    input {
        File pre
        File chrom_sizes
        RuntimeEnvironment runtime_environment
    }

    call hic.create_accessibility_track { input:
        pre = pre,
        chrom_sizes = chrom_sizes,
        runtime_environment = runtime_environment,
    }
}
