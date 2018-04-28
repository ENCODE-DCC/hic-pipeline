workflow hic {
    Array[File] fastq_files
    File restriction_sites
    File chrsz
    File reference_index

    call align { input:
        restriction = restriction_sites,
        fastqs = fastq_files,
        chrsz = chrsz,
        bwa_index = reference_index
    }

    call merge { input:
        bams = align.out_files
    }
}


task align {
    File chrsz
    File restriction
    Array[File] fastqs
    File bwa_index


    command {
        
        mkdir data && cd data && mkdir fastq && cd ..
        cd data && mkdir reference
        data_path=$(pwd)
        echo $(pwd) > data.txt
        cd ..
        cd data/fastq
        ln -s ${fastqs[0]} $(pwd)/frag_R1.fastq.gz
        ln -s ${fastqs[1]} $(pwd)/frag_R2.fastq.gz
        echo $(ls) > res2.txt
        cd ../..
        cd data/reference
        echo ${bwa_index}
        tar -xvf ${bwa_index}
        ref_path="$data_path/reference/hg38_chr19_chrM-bwa_index-GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta"
        cd ../..
        bash /opt/scripts/juicer.sh -D /opt -d $data_path -S alignonly -z $ref_path -p ${chrsz} -y ${restriction} -s MboI
    }

    output {
        Array[File] out_files = glob("data/splits/*.bam")
       
    }

    runtime {
        docker : "aidenlab/juicer:latest"
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
        docker : "aidenlab/juicer:latest"
    }
}