version 1.0

import "../../../hic.wdl" as hic

workflow test_dedup {
    input {
        File bam
        Int num_cpus
        RuntimeEnvironment runtime_environment
    }

    call hic.dedup { input:
        bam = bam,
        num_cpus = num_cpus,
        runtime_environment = runtime_environment,
    }
}
