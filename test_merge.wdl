import "hic_sub.wdl" as sub
import "hic.wdl" as hic

workflow test_merge {
    Array[File] bams
    File ref_merged
     
     
	 call sub.merge as test { input:
	 bam = bams
	 }
       
	 call hic.strip_headers as head { input:
	 bam = test.out_file
	 }  
	
	 call hic.strip_headers as ref_head { input:
	 bam = ref_merged
	 }  

	 File result = head.no_header
     File ref = ref_head.no_header
	 
	 call hic.compare_md5sum { input :
		labels = [
			'merged_bams'
		],
       
		files = [
			result			
		],

		ref_files = [
			ref	
		],
	}

}