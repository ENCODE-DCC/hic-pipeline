version 1.0

import "../../../hic.wdl" as hic

workflow test_exit_early {
    input {
        String message
        RuntimeEnvironment runtime_environment
    }

    call hic.exit_early { input:
        message = message,
        runtime_environment = runtime_environment,
    }
}
