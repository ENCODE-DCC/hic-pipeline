version 1.0

import "../../../hic.wdl" as hic

workflow test_align {
    input {
        FastqPair fastq_pair
        File idx_tar
        String ligation_site
        Int num_cpus
    }

    call hic.align { input:
        fastq_pair = fastq_pair,
        idx_tar = idx_tar,
        ligation_site = ligation_site,
        num_cpus = num_cpus,
    }
}
