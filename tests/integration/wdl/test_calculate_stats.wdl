version 1.0

import "../../../hic.wdl" as hic

workflow test_calculate_stats {
    input {
        Array[File] alignment_stats
        File bam
        File pre
        File? restriction_sites
        File chrom_sizes
        String ligation_site
        Int quality
        RuntimeEnvironment runtime_environment
    }

    call hic.calculate_stats { input:
        alignment_stats = alignment_stats,
        bam = bam,
        pre = pre,
        restriction_sites = restriction_sites,
        chrom_sizes = chrom_sizes,
        ligation_site = ligation_site,
        quality = quality,
        runtime_environment = runtime_environment,
    }
}
