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

    call merge_pairs_file { input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else process_library.library_dedup
    }


    call create_hic { input:
        pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
        #pairs_file = if defined(input_pairs) then input_pairs else dedup.out_file,
        chrsz_ = chrsz
    }

    call arrowhead { input:
        hic_file = if defined(input_hic) then input_hic else create_hic.inter_30
    }

    call hiccups{ input:
        hic_file = if defined(input_hic) then input_hic else create_hic.inter_30
    }

    output{
        # Align task outputs
        # Array[Array[File]] out_collisions = hic_sub.out_collisions
        # Array[Array[File]] out_collisions_low = hic_sub.out_collisions_low
        # Array[Array[File]] out_unmapped = hic_sub.out_unmapped
        # Array[Array[File]] out_mapq0 = hic_sub.out_mapq0
        # Array[Array[File]] out_alignable = hic_sub.out_alignable
        # Array[Array[File]] out_sort_file = hic_sub.out_sort_file
        
        # # #Merge task outputs
        # Array[File] out_merged_collisions = hic_sub.out_merged_collisions
        # Array[File] out_merged_collisions_low = hic_sub.out_merged_collisions_low
        # Array[File] out_merged_unmapped = hic_sub.out_merged_unmapped
        # Array[File] out_merged_mapq0 = hic_sub.out_merged_mapq0
        # Array[File] out_merged_align = hic_sub.out_merged_align
        
        # # #Merge sort outputs
        # Array[File] out_merge_sort = hic_sub.out_merge_sort
        # # #Dedup outputs
        # Array[File] out_dedup = hic_sub.out_dedup
        
        # # #Merge pairs file outputs
        # File out_merged_pairs = merge_pairs_file.out_file
        # # #Create hic outputs
        # File out_hic = create_hic.out_file
        # #TADs output
        # # File out_tads = tads.out_file
        # HiCCUps output
        # File out_hiccups = hiccups.out_file

        #QC outputs
        #Array[File] out_align_qc = hic_sub.out_align_qc
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
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}


task create_hic {
    File pairs_file
    File chrsz_

    command {
        /opt/scripts/common/juicer_tools pre -s inter_30.txt -g inter_30_hists.m -q 30 ${pairs_file} inter_30.hic ${chrsz_}
    }

    output {
        #/opt/scripts/common/juicer_tools pre -s inter.txt -g inter_hists.m -q 1 ${pairs_file} inter.hic ${chrsz_}
        # File inter_hic = glob('inter.hic')[0]
        File inter_30= glob('inter_30.hic')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
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
