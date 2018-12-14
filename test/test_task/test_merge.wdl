##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "../../hic.wdl" as hic

workflow test_merge {
    Array[Array[File]] bams
     
	call hic.merge as test_merge_task { input:
		bam_files = bams
    }
}