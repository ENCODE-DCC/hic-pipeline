version 1.0

import "../../../hic.wdl" as hic

workflow test_get_ligation_site_regex {
    input {
        Array[String] restriction_enzymes
    }

    call hic.get_ligation_site_regex as get_ligation_site_regex { input:
        restriction_enzymes = restriction_enzymes,
    }
}
