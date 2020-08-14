version 1.0

import "../../../hic.wdl" as hic

workflow test_exit_early {
    input {
        String message
    }

    call hic.exit_early { input:
        message = message,
    }
}
