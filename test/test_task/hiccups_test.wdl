
workflow hiccups_test {
    File input_hic
    call hiccups { input:
        hic_file = input_hic
    }
    output {
        File hiccups_out = hiccups.out_file
    } 
}

task hiccups {
    File hic_file
    command {
        java -jar /opt/scripts/common/juicer_tools.jar hiccups --ignore_sparsity ${hic_file} loops
    }
    output {
        File out_file = glob("loops/*.bedpe")[0]
    }
    runtime {
        docker: "quay.io/anacismaru/nvidia_juicer:test"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: ["us-east1-b"]
    }
}
