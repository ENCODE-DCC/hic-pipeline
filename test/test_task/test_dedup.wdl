import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_dedup {
    File merged_sort
   
    call hic.dedup as test_dedup_task { input:
    	merged_sort = merged_sort
    }
   
    output{
        File deduped = test_dedup_task.out_file
    }
}