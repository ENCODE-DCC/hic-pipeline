version 1.0

import "../../hic.wdl" as hic

workflow test_hic {
    input {
        Array[Array[Array[File]]] fastq
        String assembly_name
        String restriction_enzyme
        File restriction_sites
        File chrsz
        File reference_index
    }

    call hic.hic { input:
        fastq = fastq,
        restriction_enzyme = restriction_enzyme,
        restriction_sites = restriction_sites,
        reference_index = reference_index,
        chrsz = chrsz,
        assembly_name = assembly_name,
        no_call_loops = true,
        no_call_tads = true,
    }

    call tail_of_pairs { input:
        pairs = select_first([select_first([hic.out_pairs])[0]])
    }

    call strip_headers { input:
        bam = select_first([hic.alignable_bam])[0]
    }

    output {
        File no_header_alignable_sam = strip_headers.no_header
        File out_pairs = tail_of_pairs.no_header
        File out_dedup = select_first([hic.out_dedup])[0]
        File no_header_hic_1 = select_first([hic.out_hic_1])
        File no_header_hic_30 = select_first([hic.out_hic_30])
        File library_complexity = select_first([hic.library_complexity_stats_json])[0]
        File stats = select_first([hic.stats])[0]
        File alignments_stats = select_first([hic.alignment_stats])[0][0]
    }
}

task tail_of_pairs {
    input {
        File pairs
    }

    command {
        sed 1,5d ${pairs} > no_header.pairs
    }
    
    output {
        File no_header = glob("no_header.pairs")[0]
    }
}

task strip_headers {
    input {
        File bam
    }

    #it messes up with compare_md5.py since all the files with stripped header are having the same name
    command {
        FILE=$(basename "${bam}" ".bam")
        samtools view -h ${bam} | samtools view - > $FILE.no_header.sam
    }
    output {
        File no_header = glob("*.no_header.sam")[0]
    }
}
