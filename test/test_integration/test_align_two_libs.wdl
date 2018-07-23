##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "hic_sub.wdl" as sub
import "test/test_utils.wdl" as utils

workflow test_align {
  
	File idx_tar 		# reference bwa index tar
	Array[Array[Array[File]]] fastq = [] 	
    File chrsz          # chromosome sizes file
    File restriction    # restriction enzyme sites in the reference genome

	Array[File] ref_bams

    #determine range of scatter
    Int lib_length = if length(fastq) > 0 then length(fastq)
	
    scatter(i in range(lib_length)){
        call sub.align as test { input:
            restriction = restriction,
            fastqs = fastqs,
            chrsz = chrsz,
            idx_tar = idx_tar
        } 
    }
	Array[File] bams = test.out_file

	scatter(bam in bams){
		call utils.strip_headers as head { input:
			bam = bam
		} 
	}

	scatter(bam in ref_bams){
		call utils.strip_headers as ref_head { input:
			bam = bam
		} 
	}

	Array[File] strip_sams = head.no_header
	Array[File] ref_strip_sams = ref_head.no_header
	

    call utils.compare_md5sum { input :
		labels = [
			'align_collisions',
			'align_collisions_low_mapq',
			'align_unmapped',
			'align_mapq0',
			'align_alignable'
		],

		files = [
			strip_sams[0],
            strip_sams[1],
            strip_sams[2],
            strip_sams[3],
            strip_sams[4]
		],
    
		ref_files = [
			ref_strip_sams[0],
			ref_strip_sams[1],
			ref_strip_sams[2],
			ref_strip_sams[3],
			ref_strip_sams[4]
		],
	}
    
}