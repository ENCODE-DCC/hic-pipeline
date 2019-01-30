import "../../workflow/sub_workflow/process_library.wdl" as hic

workflow test_merge_sort {
    Array[File] sort_files
   
    call hic.merge_sort as test_merge_sort_task { input:
        sort_files_ = sort_files
    }

    output {
        File merged = test_merge_sort_task.out_file
    }
}