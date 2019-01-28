##Encode DCC Hi-C pipeline
##Author: Ana Cismaru(anacismaru@gmail.com)

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
    #scatter(i in range(lib_length)){
    
    Array[Array[File]] sub_fastq = if length(fastq) > 0 then fastq[0] else []
    Array[Array[File]] sub_input_bams = if length(input_bams) > 0 then input_bams[0] else []
    Array[File] sub_input_sort_files = if length(input_sort_files) > 0 then input_sort_files[0] else []
    File? sub_input_merged_sort = if length(input_merged_sort)>0 then input_merged_sort[0] else sub_ms



    Int fastqs_len = length(sub_fastq)
    scatter(j in range(fastqs_len)){
        call align { input:
            restriction = restriction_sites,
            fastqs = sub_fastq[j],
            chrsz = chrsz,
            idx_tar = reference_index,
            cpu = cpu 
        }
    }

    Array[File] collisions = if length(sub_input_bams)>0 then sub_input_bams[0] else align.collisions  #for separate user entry point
    Array[File] collisions_low = if length(sub_input_bams)>0 then sub_input_bams[1] else align.collisions_low_mapq
    Array[File] unmapped = if length(sub_input_bams)>0 then sub_input_bams[2] else align.unmapped
    Array[File] mapq0 = if length(sub_input_bams)>0 then sub_input_bams[3] else align.mapq0
    Array[File] alignable = if length(sub_input_bams)>0 then sub_input_bams[4] else align.alignable       
    
    Array[Array[File]] bams_to_merge = [collisions, collisions_low, unmapped, mapq0, alignable]

    scatter(i in range(5)){
        call merge { input:
            bam_files = bams_to_merge[i]
        }
    }

    # convert alignable bam to pairs to be consistent with 4DN
    call bam2pairs { input:
        bam_file = merge.merged_output[4],
        chrsz_  = chrsz
    }  
  
    call merge_sort { input:
        sort_files_ = if length(sub_input_sort_files)>0 then sub_input_sort_files else align.sort_file 
    } 

    # call align_qc { input:
    #     norm_res = align.norm_res
    # }   

    # we can collect the alignable.bam using the array merge.out_file
    call dedup { input:
        merged_sort = if defined(sub_input_merged_sort) then sub_input_merged_sort else merge_sort.out_file
    }

    #}
  
    # should be used only for multiple libraries - which we are dumping for now due to scatter within scatter limitation
    #call merge_pairs_file{ input:
    #    not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else dedup.out_file
    #}


    call create_hic { input:
        #pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
        pairs_file = if defined(input_pairs) then input_pairs else dedup.out_file,
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

task align {
    File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [read_end_id]
    File chrsz          # chromosome sizes file
    File restriction    # restriction enzyme sites in the reference genome

    Int? cpu
  

    command {     
        echo "Starting align"  
        mkdir reference
        cd reference && tar -xvf ${idx_tar}
        index_folder=$(ls)
        reference_fasta=$(ls | head -1) 
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ..
        
        usegzip=1
        curr_ostem="result"
        #HindIII site
        ligation="GATCGATC"
        file1=${fastqs[0]}
        file2=${fastqs[1]}
        #count ligations
        source /opt/scripts/common/countligations.sh
        # Align reads
        echo "Running bwa command"
        bwa mem -SP5M -t ${select_first([cpu,32])} $reference_index_path ${fastqs[0]} ${fastqs[1]} | awk -v "fname"=result -f /opt/scripts/common/chimeric_blacklist.awk
 
        # if any normal reads were written, find what fragment they correspond
        # to and store that
        echo "Running fragment"
        /opt/scripts/common/fragment.pl result_norm.txt result_frag.txt ${restriction}   
        echo $(ls)

        # qc for alignment portion
        cat *.res.txt | awk -f /opt/scripts/common/stats_sub.awk >> alignment_stats.txt

        # convert sams to bams and delete the sams
        echo "Converting sam to bam"
        samtools view -hb result_collisions.sam > collisions.bam
        rm result_collisions.sam
        samtools view -hb result_collisions_low_mapq.sam > collisions_low_mapq.bam
        rm result_collisions_low_mapq.sam
        samtools view -hb result_unmapped.sam > unmapped.bam
        rm result_unmapped.sam
        samtools view -hb result_mapq0.sam > mapq0.bam
        rm result_mapq0.sam
        samtools view -hb result_alignable.sam > alignable.bam
        rm result_alignable.sam
        #removed all sam files
        ##restriction used to be site_file
    
        # sort by chromosome, fragment, strand, and position
        sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 90% result_frag.txt > sort.txt
        if [ $? -ne 0 ]
        then
            echo "***! Failure during sort"
            exit 1x
        fi

    }

    output {
        File collisions = glob("collisions.bam")[0]
        File collisions_low_mapq = glob("collisions_low_mapq.bam")[0]
        File unmapped = glob("unmapped.bam")[0]
        File mapq0 = glob("mapq0.bam")[0]
        File alignable = glob("alignable.bam")[0]
        File sort_file = glob("sort.txt")[0]
        File norm_res = glob("result_norm.txt.res.txt")[0]
        File stats_sub_result = glob("alignment_stats.txt")[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:PIP-419-import-wdl_122a7f4a-893f-42c8-9075-7d0d256f6db0"
        cpu : "32"
        memory: "64 GB"
        disks: "local-disk 1000 HDD"
    }
}


task merge {
    Array[File] bam_files
    
    command <<<
        samtools merge merged_bam_files.bam ${sep=' ' bam_files} 
    >>>

    output {
        File merged_output= glob('merged_bam_files.bam')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "8"
        disks: "local-disk 1000 HDD"
        
    }
}

task merge_sort {
   Array[File] sort_files_

    command {
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 90% ${sep=' ' sort_files_}  > merged_sort.txt
    }

    output {
        File out_file = glob('merged_sort.txt')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}

task dedup {
    File merged_sort

    command {
        touch dups.txt
        touch optdups.txt
        touch merged_nodups.txt
        awk -f /opt/scripts/common/dups.awk ${merged_sort}
    }

    output {
        File out_file = glob('merged_nodups.txt')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "8"
        disks: "local-disk 1000 HDD"
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
        cpu : "8"
        disks: "local-disk 1000 HDD"
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
        docker : "quay.io/encode-dcc/hic-pipeline:template"
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
        docker: "quay.io/anacismaru/nvidia_juicer:test"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: ["us-east1-b"]
    }     
}
