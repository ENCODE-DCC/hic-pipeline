import "hic_sub.wdl" as sub
import "hic.wdl" as hic

workflow test_merge {
    Array[Array[File]] bams
    Array[File] ref_merged
     
     
	 call sub.merge as test { input:
	 collisions = bams[0],
	 collisions_low = bams[1],
	 unmapped = bams[2],
	 mapq0 = bams[3],
	 alignable = bams[4]
	 }

	 File collisions = test.m_collisions
	 File collisions_low = test.m_coll_low
	 File unmapped = test.m_unmap
	 File mapq0 = test.m_map
	 File aligned = test.m_align

	 call hic.strip_headers as coll { input:
		bam = collisions
	 }
	 call hic.strip_headers as low { input:
		bam = collisions_low
	 }  
	 call hic.strip_headers as unmap { input:
		bam = unmapped
	 }  
	 call hic.strip_headers as map { input:
		bam = mapq0
	 }  
	 call hic.strip_headers as align { input:
		bam = aligned
	 }    
	
	scatter(ref in ref_merged){
	call hic.strip_headers as ref_head { input:
	 		bam = ref
	 } 
	}
	  

	 Array[File] result = [coll.no_header, low.no_header, unmap.no_header, map.no_header, align.no_header]
     Array[File] ref = ref_head.no_header
	 
	 call hic.compare_md5sum { input :
		labels = [
			'merged_bams'
		],
       
		files = result,

		ref_files = ref
	}

}