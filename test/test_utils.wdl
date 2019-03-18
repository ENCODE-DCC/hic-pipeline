task strip_headers {
    File bam
    
    command {
        samtools view -h ${bam} > header.sam | samtools view > no_header.sam
    } 
    output {
        File no_header = glob("no_header.sam")[0]
    }
}

task strip_hic_header {
    File hic_file
    String? filename = basename(hic_file, ".hic")
    command {
        hic_file=${hic_file}
        matrix_start=$(python3 /opt/straw/python/get_matrix_start.py $hic_file)
        matrix_start=$((matrix_start + 1))
        # tail -c +$matrix_start $hic_file > ${filename}_no_header.hic
        tail -c 10000 $hic_file > ${filename}_no_header.hic
    }
    output {
        File no_header = glob("*_no_header.hic")[0]
    }
}
