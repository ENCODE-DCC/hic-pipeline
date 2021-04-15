version 1.0

import "./hic.wdl"

workflow hic_from_pre {
    meta {
        version: "0.1.0"
        caper_docker: "gcr.io/hic-pipeline/new-jar:latest"
        caper_singularity: "docker://encodedcc/hic-pipeline:0.1.0"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json"
    }

    input {
        File stats
        File stats_hists
        File pre
        File pre_index
        String assembly_name
        Int quality
    }

    call hic.create_hic { input:
        pre = pre,
        pre_index = pre_index,
        quality = quality,
        stats = stats,
        stats_hists = stats_hists,
        assembly_name = assembly_name,
    }

    call hic.arrowhead { input:
        hic_file = create_hic.output_hic
    }
    call hic.hiccups { input:
        hic_file = create_hic.output_hic
    }
}
