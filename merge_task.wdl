workflow hic {
    Array[File] bam_files

    call merge { input:
        bams = bam_files
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