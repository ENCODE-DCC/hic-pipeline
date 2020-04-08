version 1.0

import "../../hic.wdl" as hic

workflow test_merge_stats {
    input {
        Array[File] alignment_stats
        Array[File] library_stats
    }

    call hic.merge_stats as test_merge_stats_task { input:
        alignment_stats = alignment_stats,
        library_stats = library_stats
    }

    output {
        File merged_stats = test_merge_stats_task.merged_stats
        File merged_stats_json = test_merge_stats_task.merged_stats_json
    }
}
