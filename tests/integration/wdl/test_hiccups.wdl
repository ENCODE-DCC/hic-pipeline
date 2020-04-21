version 1.0

import "../../../hic.wdl" as hic

workflow test_hiccups {
    input {
        File input_hic
    }

    call hic.hiccups as test_hiccups_task { input:
        hic_file = input_hic
    }
    output {
        File hiccups_out = test_hiccups_task.out_file
    }
}
