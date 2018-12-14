##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
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
}