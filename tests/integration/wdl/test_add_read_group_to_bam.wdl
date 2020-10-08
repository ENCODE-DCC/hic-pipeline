version 1.0

import "../../../hic.wdl" as hic

workflow test_add_read_group_to_bam {
    input {
        File bam
        String read_group
        Int num_cpus
    }

    call hic.add_read_group_to_bam { input:
        bam = bam,
        read_group = read_group,
        num_cpus = num_cpus,
    }
}
