##Encode DCC Hi-C pipeline
import "../../workflow/sub_workflow/process_library.wdl" as sub

workflow hic {
    #User inputs 
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[Array[File]]] input_bams = [] #[lib_id[[collisions1,collisions2],[collisions_low],[unmapped],[mapq0],[alignable]], 
    Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
    Array[File] input_merged_sort = []
    Array[File] input_dedup_pairs = []
    File? input_pairs
    File? input_hic
    File? sub_ms

    String restriction_enzyme
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
            restriction_enzyme = restriction_enzyme,
            cpu = cpu,
            restriction_sites = restriction_sites
        }
    }

    call merge_pairs_file { input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else process_library.library_dedup
    }

    Array[String] qualities = ["1", "30"]
    scatter(i in range(length(qualities))) {
        call create_hic { input:
            pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
            stats = process_library.stats[0],
            stats_hists = process_library.stats_hists[0],
            chrsz_ = chrsz,
            quality = qualities[i]
        }
    }

    call arrowhead { input:
        hic_file = if defined(input_hic) then input_hic else create_hic.inter[1]
    }

    call hiccups{ input:
        hic_file = if defined(input_hic) then input_hic else create_hic.inter[1]
    }

    output{
        # Sub-workflow processing a library outputs
        Array[File] out_merged_align = process_library.alignable_bam
        Array[File] out_pairs = process_library.pairs_file
        Array[File] out_dedup = process_library.library_dedup

        # Create hic outputs
        File out_hic_1 = create_hic.inter[0]
        File out_hic_30 = create_hic.inter[1]
        
        # TADs output
        File out_tads = arrowhead.out_file
        # HiCCUps output
        File out_hiccups = hiccups.out_file

        #QC outputs
        Array[File] library_complexity = process_library.library_stats_json
        Array[File] stats = process_library.stats_json
        Array[Array[File]] alignments_stats = process_library.alignments_stats
    }

}


task merge_pairs_file{
    Array[File] not_merged_pe
    # Need to add following line back under sort command
    # ${juiceDir}/scripts/common/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt
    
    command {
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' not_merged_pe}  > merged_pairs.txt
    }
    
    output {
        File out_file = glob('merged_pairs.txt')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}


task create_hic {
    File pairs_file
    File stats
    File stats_hists
    File chrsz_
    String quality

    command {
        /opt/scripts/common/juicer_tools pre -s ${stats} -g ${stats_hists} -q ${quality} ${pairs_file} inter_${quality}.hic ${chrsz_}
    }

    output {
        File inter = glob('inter*.hic')[0]
    }

    runtime {
        cpu : "1"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
    }
}


task arrowhead {
    File hic_file

    command {
        /opt/scripts/common/juicer_tools arrowhead ${hic_file} contact_domains
    }

    output {
        File out_file = glob('contact_domains/*.bedpe')[0]
    }

    runtime {
    }
}

task hiccups{
    File hic_file
    
    command {
        java -jar /opt/scripts/common/juicer_tools.jar hiccups --ignore_sparsity ${hic_file} loops
    }
    
    output {
        File out_file = glob("loops/*.bedpe")[0]
    }
    
    runtime {
        docker: "quay.io/encode-dcc/hiccups:stripped_down"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: ["us-east1-b"]
    }     
}
