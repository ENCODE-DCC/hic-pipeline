version 1.0

import "../../../hic.wdl" as hic

workflow test_merge {
    input {
        Array[File] bams
    }

    call hic.merge { input:
        bams = bams
    }
}
