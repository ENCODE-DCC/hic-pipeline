import "../../hic.wdl" as hic

workflow test_align {
  
	File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [read_end_id]
    File chrsz          # chromosome sizes file
    File restriction    # restriction enzyme sites in the reference genome

	call hic.align as test_align_task { input:
	 	restriction = restriction,
		fastqs = fastqs,
		chrsz = chrsz,
		idx_tar = idx_tar
	}
    
    Array[File] bams = [
        test_align_task.collisions,
        test_align_task.collisions_low_mapq,
        test_align_task.unmapped,
        test_align_task.mapq0,
        test_align_task.alignable
    ]

    Int bams_len = length(bams)
    scatter(i in range(bams_len)){
        call strip_headers { input:
            bam = bams[i]
        }
    }

    output {
        File collisions = strip_headers.no_header[0]
        File collisions_low_mapq = strip_headers.no_header[1]
        File unmapped = strip_headers.no_header[2]
        File mapq0 = strip_headers.no_header[3]
        File alignable = strip_headers.no_header[4]
        
        File sort_file = test_align_task.sort_file
        File norm_res = test_align_task.norm_res
        File stats_sub_result = test_align_task.stats_sub_result
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