version 1.0

workflow delta {
    meta {
        version: "1.3.0"
        caper_docker: "encodedcc/hic-pipeline:1.3.0_delta"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.3.0_delta"
    }

    input {
        File hic
        # Should be [5000, 10000] for in-situ, [1000, 5000] otherwise
        Array[Int] resolutions = [1000, 5000, 10000]
    }

    call deploy_delta { input:
        hic = hic,
        resolutions = resolutions,
    }
}


task deploy_delta {
    input {
        File hic
        Array[Int] resolutions
        Float threshold = 0.85
        String normalization = "SCALE"
    }

    command {
        set -euo pipefail
        python \
            "$(which Deploy.py)" \
            ~{hic} \
            /opt/deploy-delta/beta-models \
            . \
            predicted \
            ~{sep="," resolutions} \
            ~{normalization} \
            ~{threshold}
        gzip -n ./*.bedpe
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
        docker: "encodedcc/hic-pipeline:1.3.0_delta"
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
