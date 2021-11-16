version 1.0

workflow delta {
    meta {
        version: "1.4.0"
        caper_docker: "encodedcc/hic-pipeline:1.4.0_delta"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.4.0_delta"
    }

    input {
        File hic
        # Should be [5000, 10000] for in-situ, [1000, 5000] otherwise
        Array[Int] resolutions = [1000, 5000, 10000]
        String docker = "encodedcc/hic-pipeline:1.4.0_delta"
    }

    scatter(resolution in resolutions) {
        call deploy_delta { input:
            hic = hic,
            resolution = resolution,
            docker = docker,
        }
    }
}


task deploy_delta {
    input {
        File hic
        Int resolution
        Float threshold = 0.85
        String normalization = "SCALE"
        String docker
    }

    command {
        set -euo pipefail
        python \
            "$(which Deploy.py)" \
            ~{hic} \
            /opt/deploy-delta/beta-models \
            . \
            predicted \
            ~{resolution} \
            ~{normalization} \
            ~{threshold}
        gzip -n ./*.bedpe
    }

    output {
        File loops = "predicted_loops_~{resolution}.bedpe.gz"
        File domains = "predicted_domains_~{resolution}.bedpe.gz"
        File stripes = "predicted_stripes_~{resolution}.bedpe.gz"
        File loop_domains = "predicted_loop_domains_~{resolution}.bedpe.gz"
    }

    runtime {
        cpu : "2"
        bootDiskSizeGb: "20"
        disks: "local-disk 100 SSD"
        docker: "~{docker}"
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
