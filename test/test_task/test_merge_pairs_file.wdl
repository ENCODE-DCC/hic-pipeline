##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "hic.wdl" as hic
import "test/test_utils.wdl" as utils

workflow test_merge_pairs_file {
    Array[File] not_merged
    File ref_merged_pairs
     
   
    call hic.merge_pairs_file as test{ input:
     not_merged_pe = not_merged
    }
        
	File result = test.out_file
     
	call utils.compare_md5sum { input :
		labels = [
			'merged_pairs'
		],
    
		files = [
			result
		],
        
		ref_files = [
			ref_merged_pairs
		],
    }

}