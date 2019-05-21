import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_create_hic {
    Array[String] ligation_junctions
    File pairs_file
    File restriction_sites

    call hic.create_hic as test_create_hic_task { input:
        pairs_file = pairs_file,
        ligation_junctions = ligation_junctions,
        restriction_sites = restriction_sites,
        quality = "30"
    }

    output {
        File hic_file = test_create_hic_task.inter
    }
}
