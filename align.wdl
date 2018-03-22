workflow hic {
    Array[File] fastq_files
    File restriction_sites
    File chrsz
    File one
    File two
    File three
    File four
    File five
    File six

    call align { input:
        restriction = restriction_sites,
        fastqs = fastq_files,
        chrsz = chrsz,
        one = one,
        two = two,
        three = three,
        four = four,
        five = five,
        six = six
    }
}


task align {
    File chrsz
    File restriction
    Array[File] fastqs
    File one
    File two
    File three
    File four
    File five
    File six

    command {
        
        mkdir data && cd data && mkdir fastq && cd ..
        cd data
        data_path=$(pwd)
        echo $(pwd) > data.txt
        cd ..
        cd data/fastq
        ln -s ${fastqs[0]} $(pwd)/frag_R1.fastq.gz
        ln -s ${fastqs[1]} $(pwd)/frag_R2.fastq.gz
        echo $(ls) > res2.txt
        cd ../..
        echo $ref_path > res.txt
        bash /opt/scripts/juicer.sh -D /opt -d $data_path -S alignonly -z ${two} -p ${chrsz} -y ${restriction} -s MboI
    }

    output {
        Array[File] out_file = glob('/splits/*.sam')
    }

    runtime {
        docker : "aidenlab/juicer:latest"
    }
}