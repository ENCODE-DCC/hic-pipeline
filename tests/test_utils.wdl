task strip_headers {
    File bam

    command {
        samtools view -h ${bam} > header.sam | samtools view > no_header.sam
    }
    output {
        File no_header = glob("no_header.sam")[0]
    }
}
