version 1.0

import "../../../hic.wdl" as hic

workflow test_merge {
    input {
        Array[File] bams
        Int num_cpus
    }

    call hic.merge { input:
        bams = bams,
        num_cpus = num_cpus,
    }
}
