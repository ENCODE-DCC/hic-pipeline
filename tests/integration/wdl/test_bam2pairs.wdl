version 1.0

import "../../hic.wdl" as hic

workflow test_bam2pairs {
    input {
        File bam
        File chrsz
    }

    call hic.bam2pairs as test { input:
        bam_file = bam,
        chrsz_ = chrsz
    }

    output {
        File pairs_no_header = test.out_file
    }
}
