import "../../workflow/main_workflow/hic.wdl" as hic
import "../../test/test_utils.wdl" as utils

workflow test_no_bam2pairs {
    #User inputs 
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    String restriction_enzyme
    File restriction_sites
    File chrsz
    File reference_index
    Boolean no_bam2pairs
    Boolean no_call_loops
    Boolean no_call_tads

    call hic.hic { input:
        fastq = fastq,
        restriction_enzyme = restriction_enzyme,
        restriction_sites = restriction_sites,
        reference_index = reference_index,
        chrsz = chrsz,
        restriction_sites = restriction_sites,
        no_call_loops = no_call_loops,
        no_call_tads = no_call_tads,
        no_bam2pairs = no_bam2pairs
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
