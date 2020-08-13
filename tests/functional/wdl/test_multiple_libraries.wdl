version 1.0

import "../../../hic.wdl" as hic

workflow test_multiple_libraries {
    input {
        Array[Array[Array[File]]] fastq
        Array[String] restriction_enzymes
        String assembly_name
        File restriction_sites
        File chrsz
        File reference_index
    }

    call hic.hic { input:
        fastq = fastq,
        chrsz = chrsz,
        restriction_sites = restriction_sites,
        reference_index = reference_index,
        restriction_enzymes = restriction_enzymes,
        assembly_name = assembly_name,
        no_call_loops = true,
        no_call_tads = true,
    }

    output{
        File alignable_bam = select_first([hic.alignable_bam])[0]
        File out_pairs = select_first([select_first([hic.out_pairs])[0]])
        File out_dedup = select_first([hic.out_dedup])[0]
        File hic_1 = select_first([hic.out_hic_1])
        File hic_30 = select_first([hic.out_hic_30])
    }
}
