version 1.0

import "../../../hic.wdl" as hic

workflow test_hic {
    input {
        Array[Array[Array[File]]] fastq
        String assembly_name
        Array[String] restriction_enzymes
        File restriction_sites
        File chrsz
        File reference_index
    }

    call hic.hic { input:
        fastq = fastq,
        restriction_enzymes = restriction_enzymes,
        restriction_sites = restriction_sites,
        reference_index = reference_index,
        chrsz = chrsz,
        assembly_name = assembly_name,
        no_call_loops = true,
        no_call_tads = true,
    }

    output {
        File alignable_bam = select_first([hic.alignable_bam])[0]
        File out_pairs = select_first([select_first([hic.out_pairs])[0]])
        File out_dedup = select_first([hic.out_dedup])[0]
        File hic_1 = select_first([hic.out_hic_1])
        File hic_30 = select_first([hic.out_hic_30])
        File library_complexity = select_first([hic.library_complexity_stats_json])[0]
        File stats = select_first([hic.stats])[0]
        File alignments_stats = select_first([hic.alignment_stats_])[0][0]
    }
}
