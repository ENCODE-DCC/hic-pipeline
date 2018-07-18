##Encode DCC Hi-C pipeline
##Author Ana Cismaru (anacismaru@gmail.com)

workflow hic_sub{
    Array[Array[File]] sub_fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[File]] sub_input_bams = []  #[collisions, collisions_low_mapq, unmapped, mapq0, alignable]
    Array[File] sub_input_sort_files = []
    File? sub_input_merged_sort #input to dedup
    
    File sub_restriction_sites
    File sub_chrsz
    File sub_reference_index
    

    #input: fastqs
    #output: bam files
    Int fastqs_len = length(sub_fastq)
    scatter(i in range(fastqs_len)){
        call align { input:
            restriction = sub_restriction_sites,
            fastqs = sub_fastq[i],
            chrsz = sub_chrsz,
            idx_tar = sub_reference_index
        }
    }
      
    Array[File] collisions = if length(sub_input_bams)>0 then sub_input_bams[0] else align.collisions  #for separate user entry point
    Array[File] collisions_low = if length(sub_input_bams)>0 then sub_input_bams[1] else align.collisions_low_mapq
    Array[File] unmapped = if length(sub_input_bams)>0 then sub_input_bams[2] else align.unmapped
    Array[File] mapq0 = if length(sub_input_bams)>0 then sub_input_bams[3] else align.mapq0
    Array[File] alignable = if length(sub_input_bams)>0 then sub_input_bams[4] else align.alignable       
    
    #input: bam files
    #output: Array of merged bam files 
    call merge { input:
        collisions = collisions,
        collisions_low = collisions_low,
        unmapped = unmapped,
        mapq0 = mapq0,
        alignable = alignable
    }
        

    #input: sort.txt 
    #output: Array of merged sort.txt
    call merge_sort { input:
        sort_files_ = if length(sub_input_sort_files)>0 then sub_input_sort_files else align.sort_file,  
    } 

    # call align_qc { input:
    #     norm_res = align.norm_res
    # }   

    # we can collect the alignable.bam using the array merge.out_file
    call dedup { input:
        merged_sort = if defined(sub_input_merged_sort) then sub_input_merged_sort else merge_sort.out_file
    }
    
    output{
        Array[File] out_collisions = align.collisions
        Array[File] out_collisions_low = align.collisions_low_mapq
        Array[File] out_unmapped = align.unmapped
        Array[File] out_mapq0 = align.mapq0
        Array[File] out_alignable = align.alignable
        Array[File] out_sort_file = align.sort_file

        File out_merged_collisions = merge.m_collisions
        File out_merged_collisions_low = merge.m_coll_low
        File out_merged_unmapped = merge.m_unmap
        File out_merged_mapq0 = merge.m_map
        File out_merged_align = merge.m_align

        File out_merge_sort = merge_sort.out_file
        File out_dedup = dedup.out_file
        #File out_align_qc = align_qc.out_file
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
        mkdir data && cd data && mkdir reference
        data_path=$(pwd)
        cd reference && tar -xvf ${idx_tar}
        index_folder=$(ls)
        cd $index_folder
        reference_fasta=$(ls | head -1) 
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ../..
       
       # Align reads
        echo "Running bwa command"
        bwa mem -SP5M -t ${select_first([cpu,4])} $reference_index_path ${fastqs[0]} ${fastqs[1]} > result.sam
       #GOOD UNTIL HERE
        
	    # chimeric takes in $name$ext
        echo "Running chimeric script"
	    awk -v "fname"=result -f /opt/scripts/common/chimeric_blacklist.awk result.sam

        # if any normal reads were written, find what fragment they correspond
 	    # to and store that
 	    echo "Running fragment"
        echo $restriction
        /opt/scripts/common/fragment.pl result_norm.txt result_frag.txt ${restriction}   
       	
        # convert sams to bams and delete the sams
        echo "Converting sam to bam"
	    samtools view -hb result_collisions.sam > collisions.bam
        samtools view -hb result_collisions_low_mapq.sam > collisions_low_mapq.bam
        samtools view -hb result_unmapped.sam > unmapped.bam
        samtools view -hb result_mapq0.sam > mapq0.bam
        samtools view -hb result_alignable.sam > alignable.bam

        #removed all sam files
        ##restriction used to be site_file
        rm result_collisions.sam result_collisions_low_mapq.sam result_unmapped.sam result_mapq0.sam result_alignable.sam
    
        # sort by chromosome, fragment, strand, and position
	    sort -T /opt/HIC_tmp -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n result_frag.txt > sort.txt
        if [ $? -ne 0 ]
    	then
            echo "***! Failure during sort"
            exit 1x
	    fi
        
    }

    output {
        File collisions = glob("data/collisions.bam")[0]
        File collisions_low_mapq = glob("data/collisions_low_mapq.bam")[0]
        File unmapped = glob("data/unmapped.bam")[0]
        File mapq0 = glob("data/mapq0.bam")[0]
        File alignable = glob("data/alignable.bam")[0]
        File sort_file = glob("data/sort.txt")[0]
        File norm_res = glob("data/result_norm.txt.res.txt")[0]
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05232018"
        cpu : 32
        memory: "64 GB"
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
        docker : "quay.io/gabdank/juicer:encode05022018"
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
        docker : "quay.io/gabdank/juicer:encode05022018"
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
        docker : "quay.io/gabdank/juicer:encode05022018"
    }
}

task align_qc {
    Array[String] norm_res
	

    command <<<
		python <<CODE	
		import json
  count_total_reads = 0
  count_unmapped = 0
  count_reg = 0
  count_norm = 0
  count_collisions = 0
  count_lowqcollisions = 0
  count_mapq0 = 0
  with open('${write_lines(norm_res)}','r') as fp:
      norm_res = fp.read().splitlines()
  for norm in norm_res:
      with open(norm,'r') as fp:
          values = fp.read().split()
          count_total_reads += int(values[0])
          count_unmapped += int(values[1])
          count_reg += int(values[2])
          count_norm += int(values[3])
          count_collisions += int(values[4])
          count_lowqcollisions += int(values[5])
          count_mapq0 += int(values[6])
      fp.close()
  data = {"Total reads": count_total_reads, "Total Unmapped": count_unmapped, "Total regular": count_reg,"Total normal": count_norm, "Total collisions": count_collisions, "Total lowqcollisions": count_lowqcollisions,"Total mapq0": count_mapq0}

  with open('align_qc.json', 'w') as outfile:
      json.dump(data, outfile)
  outfile.close()
  CODE
    >>>
    output{
        File out_file = glob("align_qc.json")[0]
    }
    runtime{
        cpu : 1
		memory : "4000 MB"
		disks : "local-disk 50 HDD" 
    }
}