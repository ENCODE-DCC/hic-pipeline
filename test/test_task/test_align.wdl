import "../../workflow/sub_workflow/process_library.wdl" as hic

workflow test_align {
    File idx_tar 		# reference bwa index tar
    Array[File] fastqs = []	# [read_end_id]
    File chrsz          # chromosome sizes file
    File restriction    # restriction enzyme sites in the reference genome
    String restriction_enzyme

    Map[String, String] restriction_enzyme_to_site = read_map("workflow/sub_workflow/restriction_enzyme_to_site.tsv")
    String ligation_site = restriction_enzyme_to_site[restriction_enzyme]

	call hic.align as test_align_task { input:
		fastqs = fastqs,
		chrsz = chrsz,
		idx_tar = idx_tar,
        ligation_site = ligation_site
	}
    
    call hic.fragment as test_fragment { input:
        restriction = restriction,
        bam_file = test_align_task.result,
        norm_res_input = test_align_task.norm_res
    }

    Array[File] bams = [
        test_fragment.collisions,
        test_fragment.collisions_low_mapq,
        test_fragment.unmapped,
        test_fragment.mapq0,
        test_fragment.alignable
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
        
        File sort_file = test_fragment.sort_file
        File norm_res = test_fragment.norm_res
        File stats_sub_result = test_fragment.stats_sub_result
        File stats_sub_result_json = test_fragment.stats_sub_result_json
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