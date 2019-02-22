import "../../workflow/main_workflow/hic.wdl" as hic

workflow test_create_hic {
    File pairs_file
    File chrsz_
    File stats
    File stats_hists

    call hic.create_hic as test_create_hic_task { input:
        pairs_file = pairs_file,
        chrsz_ = chrsz_,
        stats = stats,
        stats_hists = stats_hists,
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
