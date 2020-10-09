version 1.0

import "../../../hic.wdl" as hic

workflow test_create_hic {
    input {
        String ligation_site
        File pairs_file
        File restriction_sites
        File chrsz
        String assembly_name
    }

    call hic.create_hic as test_create_hic_task { input:
        pairs_file = pairs_file,
        ligation_site = ligation_site,
        restriction_sites = restriction_sites,
        chrsz_ = chrsz,
        quality = "30",
        assembly_name = assembly_name,
    }

    output {
        File hic_file = test_create_hic_task.inter
        File stats = test_create_hic_task.stats
        File stats_json = test_create_hic_task.stats_json
        File stats_hists = test_create_hic_task.stats_hists
    }
}
