version 1.0

import "../../../hic.wdl" as hic

workflow test_no_bam2pairs {
    input {
        Array[Array[Array[File]]] fastq = []
        String restriction_enzyme
        String assembly_name
        File restriction_sites
        File reference_index
        File chrsz
        Boolean no_bam2pairs
        Boolean no_call_loops
        Boolean no_call_tads
    }

    call hic.hic { input:
        fastq = fastq,
        restriction_enzyme = restriction_enzyme,
        restriction_sites = restriction_sites,
        reference_index = reference_index,
        restriction_sites = restriction_sites,
        chrsz = chrsz,
        no_call_loops = no_call_loops,
        no_call_tads = no_call_tads,
        no_bam2pairs = no_bam2pairs,
        assembly_name = assembly_name
    }

    Array[File?] hic_files = [hic.out_hic_1, hic.out_hic_30]

    output {
        File? inter_1 = hic_files[0]
        File? inter_30 = hic_files[1]
    }
}
