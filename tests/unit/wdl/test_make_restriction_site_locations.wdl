version 1.0

import "../../../hic.wdl" as hic

workflow test_make_restriction_site_locations {

    input {
        File reference_fasta
        String assembly_name
        String restriction_enzyme
    }

    call hic.make_restriction_site_locations { input:
        reference_fasta = reference_fasta,
        assembly_name = assembly_name,
        restriction_enzyme = restriction_enzyme,
    }
}
