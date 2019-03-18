import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_no_call_loops {
    File input_hic
    Boolean no_call_loops

    call hic.hic { input:
        input_hic = input_hic,
        no_call_loops = no_call_loops
    }
    output {
        File? domains = hic.out_tads
    }
}