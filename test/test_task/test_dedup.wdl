import "hic_sub.wdl" as sub
import "hic.wdl" as hic

workflow test_dedup {
      File merged_sort
      File ref_dedup
     
   
     call sub.dedup as test { input:
     merged_sort = merged_sort
     }
        
	File result = test.out_file
     
	 call hic.compare_md5sum { input :
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