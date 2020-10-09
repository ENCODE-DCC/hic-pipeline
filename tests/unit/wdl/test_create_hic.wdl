version 1.0

import "../../../hic.wdl" as hic

workflow test_create_hic {
    input {
        File pairs_file
        File chrsz_
        File restriction_sites
        String ligation_site
        String quality
        String assembly_name
        Array[String] normalization_methods = []
    }

    call hic.create_hic { input:
        assembly_name = assembly_name,
        chrsz_ = chrsz_,
        ligation_site = ligation_site,
        normalization_methods = normalization_methods,
        pairs_file = pairs_file,
        quality = quality,
        restriction_sites = restriction_sites,
    }
}
