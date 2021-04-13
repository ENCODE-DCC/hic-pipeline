version 1.0

import "./hic.wdl"

workflow hic_from_pre {
    input {
        File stats
        File stats_hists
        File chrsz
        File pre
        File pre_index
        String assembly_name
        Int quality
    }

    call hic.create_hic { input:
        pre = pre,
        pre_index = pre_index,
        chrsz = chrsz,
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
