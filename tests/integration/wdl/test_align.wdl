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

    output {
        File unmapped = test_fragment.unmapped
        File mapq0 = test_fragment.mapq0
        File alignable = test_fragment.alignable

        File sort_file = test_fragment.sort_file
        File norm_res = test_fragment.norm_res
        File alignment_stats = test_fragment.alignment_stats
        File alignment_stats_json = test_fragment.alignment_stats_json
    }
}
