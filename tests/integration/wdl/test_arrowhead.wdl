version 1.0

import "../../hic.wdl" as hic

workflow test_arrowhead {
    input {
        File input_hic
    }

    call hic.arrowhead as test_arrowhead_task { input:
        hic_file = input_hic
    }

    output {
        File arrowhead_out = test_arrowhead_task.out_file
    }
}
