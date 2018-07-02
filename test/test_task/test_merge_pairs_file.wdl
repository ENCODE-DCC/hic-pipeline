import "hic.wdl" as hic

workflow test_merge_pairs_file {
     Array[File] not_merged
      File ref_merged_pairs
     
   
     call hic.merge_pairs_file as test{ input:
     not_merged_pe = not_merged
     }
        
	 File result = test.out_file
     
	 call hic.compare_md5sum { input :
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