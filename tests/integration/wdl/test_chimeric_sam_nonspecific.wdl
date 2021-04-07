version 1.0

import "../../../hic.wdl" as hic

workflow test_chimeric_sam_nonspecific {
    input {
        File bam
        File ligation_count
        Int num_cpus
    }

    call hic.chimeric_sam_nonspecific { input:
        bam = bam,
        ligation_count = ligation_count,
        num_cpus = num_cpus,
    }
}
