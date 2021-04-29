version 1.0

import "../../../hic.wdl" as hic

workflow test_normalize_assembly_name {
    input {
        String assembly_name
        String normalized_assembly_name_output_path
        String assembly_is_supported_output_path
    }

    call hic.normalize_assembly_name { input:
        assembly_name = assembly_name,
        normalized_assembly_name_output_path = normalized_assembly_name_output_path,
        assembly_is_supported_output_path = assembly_is_supported_output_path,
    }
}
