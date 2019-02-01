import "../../workflow/sub_workflow/process_library.wdl" as hic

workflow test_merge {
    Array[File] bams

	call hic.merge as test_merge_task { input:
		bam_files = bams
    }
    File merged = test_merge_task.merged_output
    call strip_headers {input: bam = merged}
    output {
        File merged_no_header = strip_headers.no_header
    }
}

task strip_headers{
    File bam
    
    #it messes up with compare_md5.py since all the files with stripped header are having the same name
    command {
        FILE=$(basename "${bam}" ".bam")
        samtools view -h ${bam} | samtools view - > $FILE.no_header.sam
    }
    
    output{
        File no_header = glob("*.no_header.sam")[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
    }
}