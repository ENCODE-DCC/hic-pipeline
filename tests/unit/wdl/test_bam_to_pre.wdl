version 1.0

import "../../../hic.wdl" as hic

workflow test_bam_to_pre {
    input {
        File bam
        Int quality
        Int num_cpus
        RuntimeEnvironment runtime_environment
    }

    call hic.bam_to_pre { input:
        bam = bam,
        quality = quality,
        num_cpus = num_cpus,
        runtime_environment = runtime_environment,
    }
}
