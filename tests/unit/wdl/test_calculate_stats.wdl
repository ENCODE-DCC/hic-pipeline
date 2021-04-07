version 1.0

import "../../../hic.wdl" as hic

workflow test_calculate_stats {
    input {
        Array[File] alignment_stats
        Array[File] duplicate_counts
        File pre
        File? restriction_sites
        String ligation_site
        Int quality
    }

    call hic.calculate_stats { input:
        alignment_stats = alignment_stats,
        duplicate_counts = duplicate_counts,
        pre = pre,
        restriction_sites = restriction_sites,
        ligation_site = ligation_site,
        quality = quality,
    }
}
