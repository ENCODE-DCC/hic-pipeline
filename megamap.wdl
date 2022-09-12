version 1.0

import "./hic.wdl"

workflow megamap {
    meta {
        version: "1.15.0"
        caper_docker: "encodedcc/hic-pipeline:1.15.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.15.0"
    }

    input {
        Array[File] bigwig_files
        Array[File] hic_files
        File chrom_sizes

        # Parameters
        Int quality = 30
        Array[Int] create_hic_in_situ_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100]
        Array[Int] create_hic_intact_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]
        Boolean intact = true

        # Resource parameters
        Int? add_norm_mthreads
        Int? add_norm_save_ram
        Int? add_norm_num_cpus
        Int? add_norm_ram_gb
        Int? add_norm_disk_size_gb

        # Pipeline images
        String docker = "encodedcc/hic-pipeline:1.15.0"
        String singularity = "docker://encodedcc/hic-pipeline:1.15.0"
        String delta_docker = "encodedcc/hic-pipeline:1.15.0_delta"
        String hiccups_docker = "encodedcc/hic-pipeline:1.15.0_hiccups"
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }

    RuntimeEnvironment hiccups_runtime_environment = {
      "docker": hiccups_docker,
      "singularity": singularity
    }

    RuntimeEnvironment delta_runtime_environment = {
      "docker": delta_docker,
      "singularity": singularity
    }

    String delta_models_path = if intact then "ultimate-models" else "beta-models"
    Array[Int] delta_resolutions = if intact then [5000, 2000, 1000] else [5000, 10000]
    Array[Int] create_hic_resolutions = if intact then create_hic_intact_resolutions else create_hic_in_situ_resolutions


    call merge_bigwigs as accessibility { input:
        bigwig_files = bigwig_files,
        chrom_sizes = chrom_sizes,
        runtime_environment = runtime_environment,
    }

    call merge_stats_from_hic_files { input:
        hic_files = hic_files,
        runtime_environment = runtime_environment,
    }

    call sum_hic_files { input:
        hic_files = hic_files,
        runtime_environment = runtime_environment,
    }

    call hic.add_norm as add_norm { input:
        hic = sum_hic_files.summed_hic,
        mthreads = add_norm_mthreads,
        save_ram = add_norm_save_ram,
        quality = quality,
        juicer_tools_jar = "/opt/juicer/CPU/juicer_tools_2.17.00.jar",
        normalization_command = "addnorm2",
        num_cpus = add_norm_num_cpus,
        ram_gb = add_norm_ram_gb,
        disk_size_gb = add_norm_disk_size_gb,
        runtime_environment = runtime_environment,
    }

    call hic.arrowhead as arrowhead { input:
        hic_file = add_norm.output_hic,
        quality = quality,
        runtime_environment = runtime_environment,
    }

    if (!intact) {
        call hic.hiccups { input:
            hic_file = add_norm.output_hic,
            quality = quality,
            runtime_environment = hiccups_runtime_environment,
        }
    }

    if (intact) {
        call hic.hiccups_2 { input:
            hic = add_norm.output_hic,
            quality = quality,
            runtime_environment = hiccups_runtime_environment,
        }

        call hic.localizer as localizer_intact { input:
            hic = add_norm.output_hic,
            loops = hiccups_2.merged_loops,
            quality = quality,
            runtime_environment = runtime_environment,
        }
    }

    call hic.create_eigenvector as create_eigenvector { input:
        hic_file = add_norm.output_hic,
        chrom_sizes = chrom_sizes,
        output_filename_suffix = "_" + quality,
        runtime_environment = runtime_environment,
    }

    call hic.create_eigenvector as create_eigenvector_10kb { input:
        hic_file = add_norm.output_hic,
        chrom_sizes = chrom_sizes,
        resolution = 10000,
        output_filename_suffix = "_" + quality,
        runtime_environment = runtime_environment,
    }

    call hic.delta as delta { input:
        hic = add_norm.output_hic,
        resolutions = delta_resolutions,
        models_path = delta_models_path,
        runtime_environment = delta_runtime_environment,
    }

    call hic.localizer as localizer_delta { input:
        hic = add_norm.output_hic,
        loops = delta.loops,
        runtime_environment = runtime_environment,
    }

    call hic.slice as slice_25kb { input:
        hic_file = add_norm.output_hic,
        resolution = 25000,
        runtime_environment = runtime_environment,
    }

    call hic.slice as slice_50kb { input:
        hic_file = add_norm.output_hic,
        resolution = 50000,
        runtime_environment = runtime_environment,
    }

    call hic.slice as slice_100kb { input:
        hic_file = add_norm.output_hic,
        resolution = 100000,
        runtime_environment = runtime_environment,
    }
}


task merge_stats_from_hic_files {
    input {
        Array[File] hic_files
        Int quality = 30
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        java \
            -jar \
            /opt/merge-stats.jar \
            ~{"inter_" + quality} \
            ~{sep=" " hic_files}
        python3 \
            "$(which jsonify_stats.py)" \
            inter_~{quality}.txt \
            stats_~{quality}.json
    >>>

    output {
        File merged_stats = "inter_~{quality}.txt"
        File merged_stats_hists = "inter_~{quality}_hists.m"
        File stats_json = "stats_~{quality}.json"
    }

    runtime {
        cpu : 1
        memory: "8 GB"
        disks: "local-disk 500 HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task merge_bigwigs {
    input {
        Array[File] bigwig_files
        File chrom_sizes
        RuntimeEnvironment runtime_environment
    }

    command <<<
        bigWigMerge \
            ~{sep=" " bigwig_files} \
            combined.bedGraph
        sort -k1,1 -k2,2n combined.bedGraph > combined.sorted.bedGraph
        bedGraphToBigWig \
            combined.sorted.bedGraph \
            ~{chrom_sizes} \
            merged.bw
    >>>

    output {
        File merged_bigwig = "merged.bw"
    }

    runtime {
        cpu : 4
        memory: "32 GB"
        disks: "local-disk 500 HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task sum_hic_files {
    input {
        Array[File] hic_files
        Int num_cpus = 16
        Int ram_gb = 100
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        java \
            -jar \
            /opt/juicer/CPU/juicer_tools_2.17.00.jar \
            sum \
            --threads ~{num_cpus} \
            summed.hic \
            ~{sep=" " hic_files}
    >>>

    output {
        File summed_hic = "summed.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 500 HDD"
        memory : "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
