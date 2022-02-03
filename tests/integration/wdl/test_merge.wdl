version 1.0

import "../../../hic.wdl" as hic

workflow test_merge {
    input {
        Array[File] bams
        RuntimeEnvironment runtime_environment
    }

    call hic.merge { input:
        bams = bams
        runtime_environment = runtime_environment,
    }
}
