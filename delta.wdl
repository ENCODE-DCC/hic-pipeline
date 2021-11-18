version 1.0

workflow delta {
    meta {
        version: "1.4.0"
        caper_docker: "encodedcc/hic-pipeline:1.4.0_delta"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.4.0_delta"
    }

    input {
        File hic
        # Should be [5000, 10000] for in-situ, [1000, 5000] for intact?
        Array[Int] resolutions = [5000, 10000]
        String docker = "encodedcc/hic-pipeline:1.4.0_delta"
    }

    call deploy_delta { input:
        hic = hic,
        resolutions = resolutions,
        docker = docker,
    }

    output {
        File loops = deploy_delta.loops
        File domains = deploy_delta.domains
        File stripes = deploy_delta.stripes
        File loop_domains = deploy_delta.loop_domains
    }
}


task deploy_delta {
    input {
        File hic
        Array[Int] resolutions
        Float threshold = 0.85
        String normalization = "SCALE"
        String stem = "predicted"
        String docker
    }

    command {
        set -euo pipefail
        python \
            "$(which Deploy.py)" \
            ~{hic} \
            /opt/deploy-delta/beta-models \
            . \
            ~{stem} \
            ~{sep="," resolutions} \
            ~{normalization} \
            ~{threshold}
        gzip -n ./*.bedpe
    }

    output {
        File loops = "~{stem}_loops_merged.bedpe.gz"
        File domains = "~{stem}_domains_merged.bedpe.gz"
        File stripes = "~{stem}_stripes_merged.bedpe.gz"
        File loop_domains = "~{stem}_loop_domains_merged.bedpe.gz"
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
