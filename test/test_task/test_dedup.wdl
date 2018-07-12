##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "hic_sub.wdl" as sub
import "test/test_utils.wdl" as utils

workflow test_dedup {
    File merged_sort
    File ref_dedup
     
   
    call sub.dedup as test { input:
     merged_sort = merged_sort
    }
        
	File result = test.out_file
     
    call utils.compare_md5sum { input :
		labels = [
			'dedups'
		],
    
		files = [
			result
		],
        
		ref_files = [
			ref_dedup
		],
	}

}