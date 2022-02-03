version 1.0

import "../../../hic.wdl" as hic_

workflow test_add_norm {
    input {
        File hic
        Array[String] normalization_methods = []
        Int quality
        RuntimeEnvironment runtime_environment
    }

    call hic_.add_norm { input:
        hic = hic,
        quality = quality,
        normalization_methods = normalization_methods,
        runtime_environment = runtime_environment,
    }
}
