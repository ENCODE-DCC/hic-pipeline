version 1.0

import "../../../hic.wdl" as hic

workflow test_create_eigenvector {
    input {
        File hic_file
    }

    call hic.create_eigenvector { input:
        hic_file = hic_file,
    }
}
