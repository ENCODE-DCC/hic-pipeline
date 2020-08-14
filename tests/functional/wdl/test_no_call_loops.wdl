version 1.0

import "../../../hic.wdl" as hic

workflow test_no_call_loops {
    input {
        File input_hic
        Boolean no_call_loops
        # Added in PIP-622, not actually used by test
        Array[String] restriction_enzymes = ["MboI"]
    }

    call hic.hic { input:
        input_hic = input_hic,
        no_call_loops = no_call_loops,
        restriction_enzymes = restriction_enzymes
    }
    output {
        File? domains = hic.out_tads
    }
}
