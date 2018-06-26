import "hic_sub.wdl" as sub
import "hic.wdl" as hic

workflow test_merge_sort {
     Array[File] sort_file
      File ref_merged_sort
     
   
     call sub.merge_sort as test { input:
     sort_files_ = sort_file
     }
        
	File result = test.out_file
     call hic.compare_md5sum { input :
		labels = [
			'merged_sort'
		],
    
		files = [
			result
		],
        #TODO FIND REPLACEMENTS FOR THESE
		ref_files = [
			ref_merged_sort
		],
	}

}