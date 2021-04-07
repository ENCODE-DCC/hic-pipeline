version 1.0

import "../../../hic.wdl" as hic

workflow test_dedup {
    input {
        File bam
        Int num_cpus
    }

    call hic.dedup { input:
        bam = bam,
        num_cpus = num_cpus,
    }
}
