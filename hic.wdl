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
    
    call merge { input:
        collisions = collisions,
        collisions_low = collisions_low,
        unmapped = unmapped,
        mapq0 = mapq0,
        alignable = alignable
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
        echo "index folder ->" 
        echo $index_folder
        #cd $index_folder
        reference_fasta=$(ls | head -1) 
        echo "reference fasta ->" 
        echo $reference_fasta

        reference_folder=$(pwd)
        echo "reference_folder ->" 
        echo $reference_folder
        reference_index_path=$reference_folder/$reference_fasta
        echo "reference_index_path ->" 
        echo $reference_index_path
        cd ../..
        
        # Align reads
        echo "Running bwa command"
        bwa mem -SP5M -t ${select_first([cpu,32])} $reference_index_path ${fastqs[0]} ${fastqs[1]} | awk -v "fname"=result -f /opt/scripts/common/chimeric_blacklist.awk
        
        # chimeric takes in $name$ext
        echo "Running chimeric script"
	    

        # if any normal reads were written, find what fragment they correspond
        # to and store that
        echo "Running fragment"
        /opt/scripts/common/fragment.pl result_norm.txt result_frag.txt ${restriction}   
        echo $(ls)

        # convert sams to bams and delete the sams
        echo "Converting sam to bam"
        samtools view -hb result_collisions.sam > collisions.bam
        #rm result_collisions.sam
        samtools view -hb result_collisions_low_mapq.sam > collisions_low_mapq.bam
        #rm result_collisions_low_mapq.sam
        samtools view -hb result_unmapped.sam > unmapped.bam
        #rm result_unmapped.sam
        samtools view -hb result_mapq0.sam > mapq0.bam
        #rm result_mapq0.sam
        samtools view -hb result_alignable.sam > alignable.bam
        #rm result_alignable.sam
        #removed all sam files
        ##restriction used to be site_file
    
        # sort by chromosome, fragment, strand, and position
        sort -T /opt/HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n result_frag.txt > sort.txt
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
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:PIP-419-import-wdl_122a7f4a-893f-42c8-9075-7d0d256f6db0"
        cpu : "32"
        memory: "64 GB"
        disks: "local-disk 1000 HDD"
    }
}


task merge {
    Array[File] collisions
    Array[File] collisions_low
    Array[File] unmapped
    Array[File] mapq0 
    Array[File] alignable
     
    command <<<
        samtools merge merged_collisions.bam ${sep=' ' collisions} 
        samtools merge merged_collisions_lowmapq.bam ${sep=' ' collisions_low}
        samtools merge merged_unmapped.bam ${sep=' ' unmapped}
        samtools merge merged_mapq0.bam ${sep=' ' mapq0}
        samtools merge merged_alignable.bam ${sep=' ' alignable} 
    >>>

    output {
        File m_collisions= glob('merged_collisions.bam')[0]
        File m_coll_low = glob('merged_collisions_lowmapq.bam')[0]
        File m_unmap = glob('merged_unmapped.bam')[0]
        File m_map = glob('merged_mapq0.bam')[0]
        File m_align = glob('merged_alignable.bam')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "32"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
    }
}

task merge_sort {
   Array[File] sort_files_

    command {
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' sort_files_}  > merged_sort.txt
    }

    output {
        File out_file = glob('merged_sort.txt')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
        cpu : "32"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
        #> 8 processors
        #> a lot of memory
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
        cpu : "32"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
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
