##Encode DCC Hi-C pipeline
##Author Ana Cismaru (anacismaru@gmail.com)

workflow hic_sub{
    Array[Array[File]] sub_fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[File]] sub_input_bams = []  #[collisions, collisions_low_mapq, unampped, mapq0, alignable]
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
        #flatten bams here  
        Array[Array[File]] bams = flatten([align.out_file, sub_input_bams])  #for separate user entry point
        
        #input: bam files
        #output: Array of merged bam files
        scatter(bam in bams){
            call merge { input:
            bam = bam
            }
        }

        #input: sort.txt 
        #output: Array of merged sort.txt
        call merge_sort { input:
         sort_files_ = if length(sub_input_sort_files)>0 then sub_input_sort_files else align.sort_file,  
         }    

        # we can collect the alignable.bam using the array merge.out_file
        call dedup { input:
        #fix
         merged_sort = if defined(sub_input_merged_sort) then sub_input_merged_sort else merge_sort.out_file
         }
    
    output{
        File out_dedup = dedup.out_file
        File out_chrsz = sub_chrsz

    }
}

task align {
	File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [read_end_id]
    File chrsz          # chromosome sizes file
    File restriction    # restriction enzyme sites in the reference genome

    command {       
        mkdir data && cd data && mkdir fastq && mkdir reference
        data_path=$(pwd)
        cd fastq
        ln -s ${fastqs[0]} $(pwd)/frag_R1.fastq.gz
        ln -s ${fastqs[1]} $(pwd)/frag_R2.fastq.gz
        cd ../reference && tar -xvf ${idx_tar}
        index_folder=$(ls)
        cd $index_folder
        reference_fasta=$(ls | head -1) 
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ../..
        bash /opt/scripts/juicer.sh -D /opt -d $data_path -S alignonly -z $reference_index_path -p ${chrsz} -y ${restriction} -s MboI
    }

    output {
        File collisions = glob("data/splits/*_collisions.bam")[0]
        File collisions_low_mapq = glob("data/splits/*_collisions_low_mapq.bam")[0]
        File unampped = glob("data/splits/*_unmapped.bam")[0]
        File mapq0 = glob("data/splits/*_mapq0.bam")[0]
        File alignable = glob("data/splits/*_alignable.bam")[0]
        #TODO: reformat last 5 variables as an array or create tsv to mapping to those locations
        Array[File] out_file = [collisions, collisions_low_mapq, unampped, mapq0, alignable]
        File sort_file = glob("data/splits/*.sort.txt")[0]
   
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05232018"
        cpu : 32
        memory: "64G"
    }
}

# task test_align {
# 	File idx_tar 		# reference bwa index tar
# 	Array[File] fastqs 	# [read_end_id]
#     File chrsz          # chromosome sizes file
#     File restriction    # restriction enzyme sites in the reference genome

#     command {       
#         mkdir reference
#         cd ../reference && tar -xvf ${idx_tar}
#         index_folder=$(ls)
#         cd $index_folder
#         reference_fasta=$(ls | head -1) 
#         reference_folder=$(pwd)
#         reference_index_path=$reference_folder/$reference_fasta
       
#        # Align reads
#         echo "Running command $bwa_cmd mem -SP5M $threadstring $refSeq $file1 $file2 > ${curr_ostem}.sam" 
#         $bwa_cmd mem -SP5M $threadstring $reference_index_path $fastqs[0] $fastqs[1] > ${curr_ostem}.sam
#         if [ $? -ne 0 ]
#         then
#             echo "***! Alignment of $file1 $file2 failed."
#             exit 1
#         else                                                            
# 	    echo "(-:  Align of ${curr_ostem}.sam done successfully"
#         fi
       
       
       
       
       
#        mkdir data && cd data && mkdir fastq && mkdir reference
#         data_path=$(pwd)
#         cd fastq
#         ln -s ${fastqs[0]} $(pwd)/frag_R1.fastq.gz
#         ln -s ${fastqs[1]} $(pwd)/frag_R2.fastq.gz
#         cd ../reference && tar -xvf ${idx_tar}
#         index_folder=$(ls)
#         cd $index_folder
#         reference_fasta=$(ls | head -1) 
#         reference_folder=$(pwd)
#         reference_index_path=$reference_folder/$reference_fasta
#         cd ../..
#         bash /opt/scripts/juicer.sh -D /opt -d $data_path -S alignonly -z $reference_index_path -p ${chrsz} -y ${restriction} -s MboI
#     }

#     output {
#         File collisions = glob("data/splits/*_collisions.bam")[0]
#         File collisions_low_mapq = glob("data/splits/*_collisions_low_mapq.bam")[0]
#         File unampped = glob("data/splits/*_unmapped.bam")[0]
#         File mapq0 = glob("data/splits/*_mapq0.bam")[0]
#         File alignable = glob("data/splits/*_alignable.bam")[0]
#         #TODO: reformat last 5 variables as an array or create tsv to mapping to those locations
#         Array[File] out_file = [collisions, collisions_low_mapq, unampped, mapq0, alignable]
#         File sort_file = glob("data/splits/*.sort.txt")[0]
   
#     }

#     runtime {
#         docker : "quay.io/gabdank/juicer:encode05232018"
#         cpu : 32
#         memory: "64G"
#     }
# }

task merge {
   Array[File] bam

   command {
       samtools merge merged.bam ${sep=' ' bam}  
   }

  output {
       File out_file = glob('merged.bam')[0]
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