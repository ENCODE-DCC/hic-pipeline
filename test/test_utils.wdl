task strip_headers{
    File bam
    
    command{
        samtools view -h ${bam} > header.sam | samtools view > no_header.sam
    } 
    output{
        File no_header = glob("no_header.sam")[0]
    }
    runtime {
        docker : "quay.io/gabdank/juicer:encode05232018"
		cpu : 1
		memory : "4000 MB"
		disks : "local-disk 50 HDD"     		
	}
}
