##Encode DCC Hi-C pipeline
##Author: Ana Cismaru(anacismaru@gmail.com)

import "http://bit.ly/2RI1azD" as sub
workflow hic {
    #User inputs 
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[Array[File]]] input_bams = [] #[lib_id[[collisions1,collisions2],[collisions_low],[unmapped],[mapq0],[alignable]], 
    Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
    Array[File] input_merged_sort = []
    Array[File] input_dedup_pairs = []
    File? input_pairs
    File? input_hic
    
    File restriction_sites
    File chrsz
    File reference_index
    Int? cpu

    #determine range of scatter
    Int lib_length = if length(fastq) > 0 then length(fastq)
    else if length(input_bams) > 0 then length(input_bams) ##technically the number should be same for bams and sort_files
    else if length(input_sort_files) > 0 then length(input_sort_files)
    else length(input_merged_sort)

    #if statement includes all necesarry for creation of hic files
    #if(!defined(input_hic)){
        #call sub workflow to support input from multiple different libraries
        scatter(i in range(lib_length)){
            File? sub_ms #to deal with multiple entries
            
            call sub.hic_sub{ input:
                sub_fastq = if length(fastq) > 0 then fastq[i] else [],
                sub_input_bams = if length(input_bams) > 0 then input_bams[i] else [],
                sub_input_sort_files = if length(input_sort_files) > 0 then input_sort_files[i] else [],
                sub_input_merged_sort = if length(input_merged_sort)>0 then input_merged_sort[i] else sub_ms,

                sub_restriction_sites = restriction_sites,
                sub_chrsz = chrsz,
                sub_reference_index =reference_index,
                sub_cpu = cpu
            }
        }
   #}
    
    call merge_pairs_file{ input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else hic_sub.out_dedup
    }


    call create_hic { input:
        pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
        chrsz_ = chrsz     
    }
        
    #     #  call qc_report{ input:
    #     #  ligation = ligation,
    #     #  merged_nodups = merge_pairs_file.out_file,
    #     #  site_file = restriction_sites
    #     #  }
    # }
    
    # call tads { input:
    #     hic_file = if defined(input_hic) then input_hic else create_hic.inter_30
    # }

    # call hiccups{ input:
    #     hic_file = if defined(input_hic) then input_hic else create_hic.inter_30      
    # }
    
    output{
        #Align task outputs
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
    command{
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' not_merged_pe}  > merged_pairs.txt
    }
    output{
        File out_file = glob('merged_pairs.txt')[0]
    }
    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "32"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
    }
}


task create_hic {
    File pairs_file
    File chrsz_

    command {
        /opt/scripts/common/juicer_tools pre -s inter_30.txt -g inter_30_hists.m -q 30 ${pairs_file} inter_30.hic ${chrsz_}
    }

    output {
        # add inter_30 stuff
        #/opt/scripts/common/juicer_tools pre -s inter.txt -g inter_hists.m -q 1 ${pairs_file} inter.hic ${chrsz_}
        #File inter_hic = glob('inter.hic')[0]
         File inter_30= glob('inter_30.hic')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "32"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
    }
}

task bam2pairs {
    File bam_file
    File chrsz_ 

    command {
        /opt/pairix-0.3.6/util/bam2pairs/bam2pairs -c ${chrsz_} ${bam_file} pairix
    }

    output {
        File out_file = glob('pairix.*.pairs')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
    }
}


task tads {
    File hic_file

    command {
        /opt/scripts/common/juicer_tools arrowhead ${hic_file} contact_domains 
    }

    output {
        File out_file = glob('contact_domains/*.bedpe')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
    }
}

task hiccups{
    File hic_file
    command{
        DIR=$(dirname "${hic_file}")
        FILE=$(basename "${hic_file}")
        docker run --runtime=nvidia --entrypoint /opt/scripts/common/juicer_tools -v $DIR:/input quay.io/anacismaru/nvidia_juicer:test hiccups /input/$FILE loops
    }
    output{
        File out_file = glob("loops/*.bedpe")[0]
    }      
}
