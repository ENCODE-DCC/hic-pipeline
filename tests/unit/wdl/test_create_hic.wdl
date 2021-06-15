version 1.0

import "../../../hic.wdl" as hic

workflow test_create_hic {
    input {
        File pre
        File pre_index
        File stats
        File stats_hists
        File chrsz
        File? restriction_sites
        String assembly_name
        Array[String] normalization_methods = []
        Int quality
    }

    call hic.create_hic { input:
        pre = pre,
        pre_index = pre_index,
        stats = stats,
        stats_hists = stats_hists,
        chrsz = chrsz,
        restriction_sites = restriction_sites,
        quality = quality,
        assembly_name = assembly_name,
        normalization_methods = normalization_methods,
    }
}
