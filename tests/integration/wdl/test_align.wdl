version 1.0

import "../../../hic.wdl" as hic

workflow test_align {
    input {
        FastqPair fastq_pair
        File idx_tar
        String ligation_site
        RuntimeEnvironment runtime_environment
    }

    call hic.align { input:
        fastq_pair = fastq_pair,
        idx_tar = idx_tar,
        ligation_site = ligation_site,
        runtime_environment = runtime_environment,
    }
}
