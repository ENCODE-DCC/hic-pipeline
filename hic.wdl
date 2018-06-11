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

    
    output {
        File alignable = align.alignable
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
        File alignable = glob("data/splits/*_alignable.bam")[0]   
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05232018"
    }
}
