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

    call tail_of_pairs { input:
        pairs = test.out_file
    }

    output {
        File pairs_no_header = tail_of_pairs.no_header
    }
}


task tail_of_pairs {
    input {
        File pairs
    }

    command{
        sed 1,5d ${pairs} > no_header.pairs
    }
    
    output {
        File no_header = glob("no_header.pairs")[0]
    }
}
