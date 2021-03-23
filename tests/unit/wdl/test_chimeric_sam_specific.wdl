version 1.0

import "../../../hic.wdl" as hic

workflow test_chimeric_sam_specific {
    input {
        File bam
        File ligation_count
        File restriction_sites
        Int num_cpus
    }

    call hic.chimeric_sam_specific { input:
        bam = bam,
        ligation_count = ligation_count,
        restriction_sites = restriction_sites,
        num_cpus = num_cpus,
    }
}
