##Encode DCC Hi-C pipeline
import "../../workflow/main_workflow/hic.wdl" as hic
import "../../workflow/sub_workflow/process_library.wdl" as sub

workflow test_hic {
    #User inputs 
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[Array[File]]] input_bams = [] #[lib_id[[collisions1,collisions2],[collisions_low],[unmapped],[mapq0],[alignable]], 
    Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
    Array[File] input_merged_sort = []
    Array[File] input_dedup_pairs = []
    File? input_pairs
    File? input_hic
    File? sub_ms

    String ligation_site
    File restriction_sites
    File chrsz
    File reference_index
    Int? cpu

    #determine range of scatter
    Int lib_length = if length(fastq) > 0 then length(fastq)
    else if length(input_bams) > 0 then length(input_bams) ##technically the number should be same for bams and sort_files
    else if length(input_sort_files) > 0 then length(input_sort_files)
    else length(input_merged_sort)

    # scatter over libraries
    scatter(i in range(lib_length)) {
        call sub.process_library as process_library { input:
            sub_fastq = fastq[i],
            chrsz = chrsz,
            reference_index = reference_index,
            ligation_site = ligation_site,
            cpu = cpu,
            restriction_sites = restriction_sites
        }
    }

    call tail_of_pairs { input:
        pairs = process_library.pairs_file[0]
    }

    call strip_headers { input:
        bam = process_library.alignable_bam[0]
    }

    call hic.merge_pairs_file as merge_pairs_file { input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else process_library.library_dedup
    }

    call hic.create_hic as create_hic { input:
        pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
        stats = process_library.stats[0],
        stats_hists = process_library.stats_hists[0],
        chrsz_ = chrsz
    }

    call strip_header { input:
        hic_file = create_hic.inter_30
    }

    output{
        File no_header_alignable_sam = strip_headers.no_header
        File out_pairs = tail_of_pairs.no_header
        File out_dedup = process_library.library_dedup[0]
        File no_header_hic = strip_header.no_header

        #QC outputs
        File library_complexity = process_library.library_stats_json[0]
        File stats = process_library.stats_json[0]
        File alignments_stats = process_library.alignments_stats[0][0]
    }
}

task tail_of_pairs{
    File pairs

    command{
        sed 1,5d ${pairs} > no_header.pairs
    }
    
    output {
        File no_header = glob("no_header.pairs")[0]
    }
}

task strip_headers{
    File bam

    #it messes up with compare_md5.py since all the files with stripped header are having the same name
    command {
        FILE=$(basename "${bam}" ".bam")
        samtools view -h ${bam} | samtools view - > $FILE.no_header.sam
    }
    
    output {
        File no_header = glob("*.no_header.sam")[0]
    }
}

task strip_header {
    File hic_file

    command {
        hic_file=${hic_file}
        matrix_start=$(python3 /opt/straw/python/get_matrix_start.py $hic_file)
        tail -c +$matrix_start $hic_file > no_header.hic
        echo $(wc -c no_header.hic)
    }

    output {
        File no_header = glob("no_header.hic")[0]
    }
}