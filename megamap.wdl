version 1.0

import "./hic.wdl"

workflow megamap {
    meta {
        version: "1.15.1"
        caper_docker: "encodedcc/hic-pipeline:1.15.1"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.15.1"
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
        Boolean no_delta = true
        Int localizer_resolution = 100

        # Resource parameters
        Int? add_norm_mthreads
        Int? add_norm_save_ram
        Int? add_norm_num_cpus
        Int? add_norm_ram_gb
        Int? add_norm_disk_size_gb

        Int? create_eigenvector_disk_size_gb
        Int? create_eigenvector_ram_gb

        Int? delta_disk_size_gb
        Int? delta_num_gpus
        Int? delta_ram_gb

        Int? slice_disk_size_gb
        Int? slice_num_cpus
        Int? slice_ram_gb

        Int? arrowhead_disk_size_gb
        Int? arrowhead_num_cpus
        Int? arrowhead_ram_gb

        Int? hiccups_2_disk_size_gb
        Int? hiccups_2_num_gpus
        Int? hiccups_2_num_cpus
        Int? hiccups_2_ram_gb

        Int? localizer_disk_size_gb

        Int? merge_bigwigs_disk_size_gb
        Int? merge_bigwigs_num_cpus
        Int? merge_bigwigs_ram_gb

        Int? sum_hic_files_disk_size_gb
        Int? sum_hic_files_ram_gb
        Int? sum_hic_files_num_cpus
        Int? sum_hic_files_num_threads

        # Pipeline images
        String docker = "encodedcc/hic-pipeline:1.15.1"
        String singularity = "docker://encodedcc/hic-pipeline:1.15.1"
        String delta_docker = "encodedcc/hic-pipeline:1.15.1_delta"
        String hiccups_docker = "encodedcc/hic-pipeline:1.15.1_hiccups"
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
        disk_size_gb = merge_bigwigs_disk_size_gb,
        num_cpus = merge_bigwigs_num_cpus,
        ram_gb = merge_bigwigs_ram_gb,
    }

    call merge_stats_from_hic_files { input:
        hic_files = hic_files,
        runtime_environment = runtime_environment,
    }

    call sum_hic_files { input:
        hic_files = hic_files,
        disk_size_gb = sum_hic_files_disk_size_gb,
        ram_gb = sum_hic_files_ram_gb,
        num_cpus = sum_hic_files_num_cpus,
        num_threads = sum_hic_files_num_threads,
        runtime_environment = runtime_environment,
    }

    call hic.add_norm as add_norm { input:
        hic = sum_hic_files.summed_hic,
        mthreads = add_norm_mthreads,
        save_ram = add_norm_save_ram,
        quality = quality,
        juicer_tools_jar = "/opt/juicer/CPU/juicer_tools.2.20.00.jar",
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
        disk_size_gb = arrowhead_disk_size_gb,
        num_cpus = arrowhead_num_cpus,
        ram_gb = arrowhead_ram_gb,
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
            disk_size_gb = hiccups_2_disk_size_gb,
            num_cpus = hiccups_2_num_cpus,
            num_gpus = hiccups_2_num_gpus,
            ram_gb = hiccups_2_ram_gb,
        }

        call hic.localizer as localizer_intact { input:
            hic = add_norm.output_hic,
            loops = hiccups_2.merged_loops,
            localizer_resolution = localizer_resolution,
            localizer_window = 10,
            quality = quality,
            disk_size_gb = localizer_disk_size_gb,
            runtime_environment = runtime_environment,
        }
    }

    call hic.create_eigenvector as create_eigenvector { input:
        hic_file = add_norm.output_hic,
        chrom_sizes = chrom_sizes,
        output_filename_suffix = "_" + quality,
        runtime_environment = runtime_environment,
        disk_size_gb = create_eigenvector_disk_size_gb,
        ram_gb = create_eigenvector_ram_gb,
    }

    call hic.create_eigenvector as create_eigenvector_10kb { input:
        hic_file = add_norm.output_hic,
        chrom_sizes = chrom_sizes,
        resolution = 10000,
        output_filename_suffix = "_" + quality,
        runtime_environment = runtime_environment,
        disk_size_gb = create_eigenvector_disk_size_gb,
        ram_gb = create_eigenvector_ram_gb,
    }

    if (!no_delta) {
        call hic.delta as delta { input:
            hic = add_norm.output_hic,
            resolutions = delta_resolutions,
            models_path = delta_models_path,
            runtime_environment = delta_runtime_environment,
            disk_size_gb = delta_disk_size_gb,
            ram_gb = delta_ram_gb,
            gpu_count = delta_num_gpus,
        }

        call hic.localizer as localizer_delta { input:
            hic = add_norm.output_hic,
            loops = delta.loops,
            localizer_resolution = localizer_resolution,
            localizer_window = 10,
            disk_size_gb = localizer_disk_size_gb,
            runtime_environment = runtime_environment,
        }
    }

    call hic.slice as slice_25kb { input:
        hic_file = add_norm.output_hic,
        resolution = 25000,
        runtime_environment = runtime_environment,
        disk_size_gb = slice_disk_size_gb,
        num_cpus = slice_num_cpus,
        ram_gb = slice_ram_gb,
    }

    call hic.slice as slice_50kb { input:
        hic_file = add_norm.output_hic,
        resolution = 50000,
        runtime_environment = runtime_environment,
        disk_size_gb = slice_disk_size_gb,
        num_cpus = slice_num_cpus,
        ram_gb = slice_ram_gb,
    }

    call hic.slice as slice_100kb { input:
        hic_file = add_norm.output_hic,
        resolution = 100000,
        runtime_environment = runtime_environment,
        disk_size_gb = slice_disk_size_gb,
        num_cpus = slice_num_cpus,
        ram_gb = slice_ram_gb,
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
        Int disk_size_gb = 500
        Int num_cpus = 4
        Int ram_gb = 32
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
        cpu : "~{num_cpus}"
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task sum_hic_files {
    input {
        Array[File] hic_files
        Int disk_size_gb = 500
        Int num_cpus = 16
        Int num_threads = 8
        Int ram_gb = 100
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        java \
            -jar \
            /opt/juicer/CPU/juicer_tools.2.20.00.jar \
            sum \
            --threads ~{num_threads} \
            summed.hic \
            ~{sep=" " hic_files}
    >>>

    output {
        File summed_hic = "summed.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory : "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
