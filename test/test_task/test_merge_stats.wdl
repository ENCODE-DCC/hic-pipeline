#CAPER docker quay.io/encode-dcc/hic-pipeline:template

import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_merge_stats {
    Array[File] alignment_stats
    Array[File] library_stats

    call hic.merge_stats as test_merge_stats_task { input:
        alignment_stats = alignment_stats,
        library_stats = library_stats
    }

    output {
        File merged_stats = test_merge_stats_task.merged_stats
        File merged_stats_json = test_merge_stats_task.merged_stats_json
    }
}
