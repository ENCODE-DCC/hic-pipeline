version 1.0

import "../../hic.wdl" as hic

workflow test_merge {
    input {
        Array[File] bams
    }

	call hic.merge as test_merge_task { input:
		bam_files = bams
    }

    output {
        File merged_no_header = test_merge_task.merged_output
    }
}
