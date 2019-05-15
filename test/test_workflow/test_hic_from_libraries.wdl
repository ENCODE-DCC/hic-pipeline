import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_hic_from_libraries {
    Array[File] input_dedup_pairs
    Array[File] library_stats
    Array[File] library_stats_hists
    File chrsz
    File restriction_sites
    String restriction_enzyme
    Boolean no_call_loops
    Boolean no_call_tads

    call hic.hic { input:
        input_dedup_pairs = input_dedup_pairs,
        library_stats = library_stats,
        library_stats_hists = library_stats_hists,
<<<<<<< HEAD
        chrsz = chrsz,
        restriction_enzyme = restriction_enzyme,
=======
        input_ligation_junctions = ligation_junctions,
>>>>>>> fix hic_from_libraries test
        restriction_sites = restriction_sites,
        no_call_loops = no_call_loops,
        no_call_tads = no_call_tads
    }

    output {
        File? inter_1_no_header = hic.out_hic_1
        File? inter_30_no_header = hic.out_hic_30
    }
}
