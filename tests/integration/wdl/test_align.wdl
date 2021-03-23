version 1.0

import "../../../hic.wdl" as hic

workflow test_align {
    input {
        File idx_tar
        Array[File] fastqs
        String ligation_site
    }

    call hic.align { input:
        fastqs = fastqs,
        idx_tar = idx_tar,
        ligation_site = ligation_site
    }
}
