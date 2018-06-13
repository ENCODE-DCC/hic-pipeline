workflow hic {
    Array[Array[File]] fastq_files
    File restriction_sites
    File chrsz
    File reference_index


    Int fastqs_len = length(fastq_files)
    # TODO
    # we need to support scatterring over libraries as well
    # TODO
    scatter(i in range(fastqs_len)){
        call align { input:
            restriction = restriction_sites,
            fastqs = fastq_files[i],
            chrsz = chrsz,
            idx_tar = reference_index
        }
    }

    Array[Array[File]] bam_files = [
       align.collisions,
       align.collisions_low_mapq,
       align.unampped,
       align.mapq0,
       align.alignable
   ]
    Int bams_len = length(bam_files)
    
     scatter(i in range(bams_len)){
        call merge { input:
           bams = bam_files[i]
        }
    }
     
     call merge_sort { input:
       sort_files = align.sort_file,  
    }
   
    # we can collect the alignable.bam using the array merge.out_file
    call dedup { input:
       merged_sort = merge_sort.out_file
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
        File sort_file = glob("data/splits/*.sort.txt")[0]
   
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05232018"
        cpu : 32
        memory: "64G"
    }
}

task merge {
   Array[File] bams

   command {
       samtools merge merged.bam ${sep=' ' bams}  
   }

  output {
       File out_file = glob('merged.bam')[0]
   }

  runtime {
       docker : "quay.io/gabdank/juicer:encode05022018"
   }
}

task merge_sort {
   Array[File] sort_files

   command {
       sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' sort_files}  > merged_sort.txt
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
