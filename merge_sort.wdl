workflow hic {
    Array[File] sort_files

    call merge { input:
        sort_files = sort_files
    }
}

task merge {
    Array[File] sort_files

    command {
        samtools merge merged.bam ${sep=' ' bams}   
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' sort_files}  > $outputdir/merged_sort.txt

    }

    output {
        File out_file = glob('merge_sorted.txt')[0]
    }

    runtime {
        docker : "aidenlab/juicer:latest"

        > 8 processors
        > a lot of memory
    }
}