import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_create_hic {
    Array[String] ligation_junctions
    File pairs_file
    File chrsz_
    File restriction_sites

    call hic.create_hic as test_create_hic_task { input:
        pairs_file = pairs_file,
        use_juicer_chrsz = true,
        ligation_junctions = ligation_junctions,
        restriction_sites = restriction_sites,
        quality = "30"
    }

    File hic_file = test_create_hic_task.inter
    
    call strip_header { input:
        hic_file = hic_file
    }
    
    output {
        File no_header = strip_header.no_header
    }
}

task strip_header {
    File hic_file
    command {
        hic_file=${hic_file}
        matrix_start=$(python3 /opt/straw/python/get_matrix_start.py $hic_file)
        matrix_start=$((matrix_start + 1))
        # tail -c +$matrix_start $hic_file > no_header.hic
        tail -c 10000 $hic_file > no_header.hic
    }
    output {
        File no_header = glob("no_header.hic")[0]
    }
}
