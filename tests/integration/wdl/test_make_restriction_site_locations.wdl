version 1.0

import "../../../make_restriction_site_locations.wdl"

workflow test_make_restriction_site_locations {

    input {
        File reference_fasta
        String assembly_name
        String restriction_enzyme
        RuntimeEnvironment runtime_environment
    }

    call make_restriction_site_locations.make_restriction_site_locations_ { input:
        reference_fasta = reference_fasta,
        assembly_name = assembly_name,
        restriction_enzyme = restriction_enzyme,
        runtime_environment = runtime_environment,
    }
}
