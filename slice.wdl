version 1.0

import "./hic.wdl"


workflow slice {
    meta {
        version: "1.1.0"
        caper_docker: "encodedcc/hic-pipeline:1.1.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.1.0"
    }

    input {
        File hic_file
    }

    call hic.slice as run_slice { input:
        hic_file = hic_file
    }
}
