version 1.0

import "../../../hic.wdl" as hic

workflow test_arrowhead {
    input {
        File input_hic
        Boolean ignore_sparsity
    }

    call hic.arrowhead as test_arrowhead_task { input:
        hic_file = input_hic,
        ignore_sparsity = ignore_sparsity,
    }

    output {
        File arrowhead_out = test_arrowhead_task.out_file
    }
}
