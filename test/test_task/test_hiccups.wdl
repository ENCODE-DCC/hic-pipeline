import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_hiccups {
    File input_hic
    call hic.hiccups as test_hiccups_task { input:
        hic_file = input_hic
    }
    output {
        File hiccups_out = test_hiccups_task.out_file
    }
}
