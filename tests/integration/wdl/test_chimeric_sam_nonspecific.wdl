version 1.0

import "../../../hic.wdl" as hic

workflow test_chimeric_sam_nonspecific {
    input {
        File bam
        File ligation_count
        Boolean single_ended
        Int num_cpus
    }

    call hic.chimeric_sam_nonspecific { input:
        bam = bam,
        ligation_count = ligation_count,
        single_ended = single_ended,
        num_cpus = num_cpus,
    }
}
