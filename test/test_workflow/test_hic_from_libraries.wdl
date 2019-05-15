import "../../workflow/main_workflow/hic.wdl" as hic
import "../../test/test_utils.wdl" as utils

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
        chrsz = chrsz,
        restriction_enzyme = restriction_enzyme,
        restriction_sites = restriction_sites,
        no_call_loops = no_call_loops,
        no_call_tads = no_call_tads
    }

    Array[File?] hic_files = [hic.out_hic_1, hic.out_hic_30]

    scatter(hic_file in hic_files) {
        call utils.strip_hic_header as strip_hic_header { input:
            hic_file = hic_file
        }
    }

    output {
        File? inter_1_no_header = strip_hic_header.no_header[0]
        File? inter_30_no_header = strip_hic_header.no_header[1]
    }
}
