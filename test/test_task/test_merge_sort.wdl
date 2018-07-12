##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "hic_sub.wdl" as sub
import "test/test_utils.wdl" as utils

workflow test_merge_sort {
    Array[File] sort_file
    File ref_merged_sort
     
   
    call sub.merge_sort as test { input:
     sort_files_ = sort_file
    }
        
	File result = test.out_file
    call utils.compare_md5sum { input :
		labels = [
			'merged_sort'
		],
    
		files = [
			result
		],
        
		ref_files = [
			ref_merged_sort
		],
	}

}