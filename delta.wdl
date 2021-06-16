version 1.0

workflow delta {
    meta {
        version: "0.7.0"
        caper_docker: "gcr.io/hic-pipeline/delta"
        caper_singularity: "docker://gcr.io/hic-pipeline/delta"
    }

    input {
        File hic
        Array[File] models
        # in-situ or intact
        String mode = "in-situ"
        # Overrides the detault resolutions that are based on the mdoe
        Array[Int] resolutions
        String? normalization
    }

    Array[Int] resolutions_ =
        if (defined(resolutions)) then resolutions
        else if (mode == "in-situ") then [5000, 10000]
        else [1000, 5000]

    call deploy_delta { input:
        hic = hic,
        models = models,
        resolutions = resolutions_,
        normalization = normalization,
    }
}


task deploy_delta {
    input {
        File hic
        Array[File] models
        Array[Int] resolutions
        String normalization = "SCALE"
    }

    command {
        set -euo pipefail
        mkdir models
        cp ~{sep=" " models} models
        python \
            "$(which Deploy.py)" \
            ~{hic} \
            models \
            . \
            predicted \
            ~{sep="," resolutions} \
            ~{normalization}
        gzip -n loops/*.bedpe
    }

    output {
        Array[File] loops = glob("predicted_loops_*")
        Array[File] domains = glob("predicted_domains_*")
        Array[File] stripes = glob("predicted_stripes_*")
        Array[File] loop_domains = glob("predicted_loop_domains_*")
    }

    runtime {
        cpu : "2"
        bootDiskSizeGb: "20"
        disks: "local-disk 100 SSD"
        # docker: "encodedcc/hic-pipeline:0.7.0_delta"
        docker: "gcr.io/hic-pipeline/delta"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        memory: "32 GB"
        zones: [
            "us-central1-c",
            "us-central1-f",
            "us-east1-b",
            "us-east1-c",
            "us-west1-a",
            "us-west1-b",
        ]
    }
}
