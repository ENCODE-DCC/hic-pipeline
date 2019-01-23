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

	output {
        File collisions = test_align_task.collisions
        File collisions_low_mapq = test_align_task.collisions_low_mapq
        File unmapped = test_align_task.unmapped
        File mapq0 = test_align_task.mapq0
        File alignable = test_align_task.alignable
        File sort_file = test_align_task.sort_file
        File norm_res = test_align_task.norm_res
        File stats_sub_result = test_align_task.stats_sub_result
    }
}