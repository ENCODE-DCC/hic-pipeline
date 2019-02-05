import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_arrowhead {
    File input_hic
    
    call hic.arrowhead as test_arrowhead_task { input:
        hic_file = input_hic
    }
    
    output {
        File arrowhead_out = test_arrowhead_task.out_file
    } 
}
