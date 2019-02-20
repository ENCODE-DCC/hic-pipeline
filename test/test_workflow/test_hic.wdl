##Encode DCC Hi-C pipeline
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
        call sub.process_library { input:
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

    call merge_pairs_file { input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else process_library.library_dedup
    }

    call create_hic { input:
        pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
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
        File library_complexity = process_library.library_stats[0]
        File stats = process_library.stats[0]
        File alignments_stats = process_library.alignments_stats[0][0]
    }
}

task merge_pairs_file{
    Array[File] not_merged_pe

    command {
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' not_merged_pe}  > merged_pairs.txt
    }
    
    output {
        File out_file = glob('merged_pairs.txt')[0]
    }
}


task create_hic {
    File pairs_file
    File chrsz_

    command {
        /opt/scripts/common/juicer_tools pre -s inter_30.txt -g inter_30_hists.m -q 30 ${pairs_file} inter_30.hic ${chrsz_}
    }

    output {
        File inter_30= glob('inter_30.hic')[0]
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