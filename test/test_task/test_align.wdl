version 1.0

import "../../hic.wdl" as hic

workflow test_align {
    input {
        File idx_tar
        Array[File] fastqs
        File chrsz
        File restriction
        String restriction_enzyme
    }

    # Can't read map from hic.wdl, so hard coded. Corresponse to MboI
    String ligation_site = "GATCGATC"

	call hic.align as test_align_task { input:
		fastqs = fastqs,
		chrsz = chrsz,
		idx_tar = idx_tar,
        ligation_site = ligation_site
	}
    
    call hic.fragment as test_fragment { input:
        restriction = restriction,
        bam_file = test_align_task.result,
        norm_res_input = test_align_task.norm_res
    }

    Array[File] bams = [
        test_fragment.collisions,
        test_fragment.collisions_low_mapq,
        test_fragment.unmapped,
        test_fragment.mapq0,
        test_fragment.alignable
    ]

    Int bams_len = length(bams)
    scatter(i in range(bams_len)){
        call strip_headers { input:
            bam = bams[i]
        }
    }

    output {
        File collisions = strip_headers.no_header[0]
        File collisions_low_mapq = strip_headers.no_header[1]
        File unmapped = strip_headers.no_header[2]
        File mapq0 = strip_headers.no_header[3]
        File alignable = strip_headers.no_header[4]
        
        File sort_file = test_fragment.sort_file
        File norm_res = test_fragment.norm_res
        File alignment_stats = test_fragment.alignment_stats
        File alignment_stats_json = test_fragment.alignment_stats_json
    }
}

task strip_headers {
    input {
        File bam
    }

    #it messes up with compare_md5.py since all the files with stripped header are having the same name
    command {
        FILE=$(basename "${bam}" ".bam")
        samtools view -h ${bam} | samtools view - > $FILE.no_header.sam
    }
    
    output{
        File no_header = glob("*.no_header.sam")[0]
    }
}
