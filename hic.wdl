version 1.0

struct FastqPair {
    File read_1
    File? read_2
    String? read_group
}

struct BamAndLigationCount {
    File bam
    File ligation_count
    Boolean single_ended
}

struct RuntimeEnvironment {
    String docker
    String singularity
}

workflow hic {
    meta {
        version: "1.15.1"
        caper_docker: "encodedcc/hic-pipeline:1.15.1"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.15.1"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json"
        description: "ENCODE Hi-C pipeline, see https://github.com/ENCODE-DCC/hic-pipeline for details."
    }

    input {
        # Main entrypoint, need to specify all five of these values except read_groups when running from fastqs
        Array[Array[FastqPair]] fastq = []
        Array[String] restriction_enzymes = []
        String? ligation_site_regex
        File? restriction_sites
        File? chrsz
        File? reference_index

        # Entrypoint for loop and TAD calls
        File? input_hic

        # Parameters controlling delta calls
        Boolean no_delta = false

        # Parameters
        Boolean intact = false
        Array[String] normalization_methods = []
        Array[Int] create_hic_in_situ_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100]
        Array[Int] create_hic_intact_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]

        # Feature flags
        Boolean no_pairs = false
        Boolean no_call_loops = false
        Boolean no_call_tads = false
        Boolean no_eigenvectors = false
        Boolean no_slice = false

        # Resource parameters
        Int align_num_cpus = 32
        Int align_ram_gb_in_situ = 64
        Int align_ram_gb_intact = 88
        Int align_disk_size_gb_in_situ = 1000
        Int align_disk_size_gb_intact = 1500
        Int chimeric_sam_nonspecific_disk_size_gb = 7000
        Int chimeric_sam_specific_disk_size_gb = 1500
        Int dedup_ram_gb_in_situ = 32
        Int dedup_ram_gb_intact = 48
        Int dedup_disk_size_gb_in_situ = 5000
        Int dedup_disk_size_gb_intact = 7500
        Int? create_hic_num_cpus
        Int? create_hic_ram_gb
        Int? create_hic_juicer_tools_heap_size_gb
        Int? create_hic_disk_size_gb
        Int? add_norm_num_cpus
        Int? add_norm_ram_gb
        Int? add_norm_disk_size_gb
        Int? create_accessibility_track_ram_gb
        Int? create_accessibility_track_disk_size_gb
        String assembly_name = "undefined"

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

    Int align_ram_gb = if intact then align_ram_gb_intact else align_ram_gb_in_situ
    Int align_disk_size_gb = if intact then align_disk_size_gb_intact else align_disk_size_gb_in_situ
    Int dedup_ram_gb = if intact then dedup_ram_gb_intact else dedup_ram_gb_in_situ
    Int dedup_disk_size_gb = if intact then dedup_disk_size_gb_intact else dedup_disk_size_gb_in_situ
    String delta_models_path = if intact then "ultimate-models" else "beta-models"
    Array[Int] delta_resolutions = if intact then [5000, 2000, 1000] else [5000, 10000]
    Array[Int] create_hic_resolutions = if intact then create_hic_intact_resolutions else create_hic_in_situ_resolutions

    # Default MAPQ thresholds for generating .hic contact maps
    Array[Int] DEFAULT_HIC_QUALITIES = [1, 30]
    Boolean is_nonspecific = length(restriction_enzymes) > 0 && restriction_enzymes[0] == "none"

    if (!defined(input_hic)) {
        if (!defined(ligation_site_regex)) {
            call get_ligation_site_regex { input:
                restriction_enzymes = restriction_enzymes,
                runtime_environment = runtime_environment,
            }
        }

        String ligation_site = select_first([ligation_site_regex, get_ligation_site_regex.ligation_site_regex])

        if (!is_nonspecific && !defined(restriction_sites)) {
            call exit_early { input:
                message = "Must provide restriction sites file if enzyme is not `none`",
                runtime_environment = runtime_environment,
            }
        }
    }

    call normalize_assembly_name { input:
        assembly_name = assembly_name,
        runtime_environment = runtime_environment,
    }

    scatter(i in range(length(fastq))) {
        Array[FastqPair] replicate = fastq[i]
        scatter(fastq_pair in replicate) {
            call align { input:
                fastq_pair = fastq_pair,
                idx_tar = select_first([reference_index]),
                ligation_site = select_first([ligation_site]),
                num_cpus = align_num_cpus,
                ram_gb = align_ram_gb,
                disk_size_gb = align_disk_size_gb,
                runtime_environment = runtime_environment,
            }
        }

        if (is_nonspecific) {
            scatter(bam_and_ligation_count in align.bam_and_ligation_count) {
                call chimeric_sam_nonspecific { input:
                    bam = bam_and_ligation_count.bam,
                    ligation_count = bam_and_ligation_count.ligation_count,
                    single_ended = bam_and_ligation_count.single_ended,
                    disk_size_gb = chimeric_sam_nonspecific_disk_size_gb,
                    runtime_environment = runtime_environment,
                }
            }
        }

        if (!is_nonspecific) {
            scatter(bam_and_ligation_count in align.bam_and_ligation_count) {
                call chimeric_sam_specific { input:
                    bam = bam_and_ligation_count.bam,
                    ligation_count = bam_and_ligation_count.ligation_count,
                    restriction_sites = select_first([restriction_sites]),
                    single_ended = bam_and_ligation_count.single_ended,
                    disk_size_gb = chimeric_sam_specific_disk_size_gb,
                    runtime_environment = runtime_environment,
                }
            }
        }

        call merge { input:
            bams = flatten(
                select_all(
                    [
                        chimeric_sam_specific.output_bam,
                        chimeric_sam_nonspecific.output_bam,
                    ]
                )
            ),
            output_bam_filename = "merged_" + i,
            runtime_environment = runtime_environment,
        }

        call dedup { input:
            bam = merge.bam,
            ram_gb = dedup_ram_gb,
            disk_size_gb = dedup_disk_size_gb,
            runtime_environment = runtime_environment,
        }

        call bam_to_pre as bam_to_pre_for_stats { input:
            bam = dedup.deduped_bam,
            quality = 1,
            output_filename_suffix = "_lib" + i,
            runtime_environment = runtime_environment,
        }

        call calculate_stats as calculate_stats_on_library { input:
            alignment_stats = flatten(
                select_all([chimeric_sam_specific.stats, chimeric_sam_nonspecific.stats])
            ),
            bam = dedup.deduped_bam,
            pre = bam_to_pre_for_stats.pre,
            restriction_sites = restriction_sites,
            chrom_sizes = select_first([chrsz]),
            ligation_site = select_first([ligation_site]),
            output_filename_suffix = "_lib" + i,
            single_ended = align.bam_and_ligation_count[0].single_ended,
            runtime_environment = runtime_environment,
        }
    }

    if (!defined(input_hic)) {
        call merge as merge_replicates { input:
            bams = dedup.deduped_bam,
            runtime_environment = runtime_environment,
        }
        # convert alignable bam to pairs to be consistent with 4DN
        if ( !no_pairs && defined(chrsz)) {
            call bam_to_pre as bam_to_pre_mapq0 { input:
                bam = merge_replicates.bam,
                quality = 0,
                runtime_environment = runtime_environment,
            }
            call pre_to_pairs { input:
                pre = bam_to_pre_mapq0.pre,
                chrom_sizes  = select_first([chrsz]),
                runtime_environment = runtime_environment,
            }
        }
    }

    Array[Int] qualities = if !defined(input_hic) then DEFAULT_HIC_QUALITIES else []
    scatter(i in range(length(qualities))) {
        call bam_to_pre { input:
            bam = select_first([merge_replicates.bam]),
            quality = qualities[i],
            runtime_environment = runtime_environment,
        }

        if (intact && qualities[i] == 30) {
            call create_accessibility_track { input:
                pre = bam_to_pre.pre,
                chrom_sizes = select_first([chrsz]),
                ram_gb = create_accessibility_track_ram_gb,
                disk_size_gb = create_accessibility_track_disk_size_gb,
                runtime_environment = runtime_environment,
            }
        }

        call calculate_stats { input:
            alignment_stats = flatten(
                select_all(
                    flatten(
                        [chimeric_sam_specific.stats, chimeric_sam_nonspecific.stats]
                    )
                )
            ),
            bam = select_first([merge_replicates.bam]),
            pre = bam_to_pre.pre,
            restriction_sites = restriction_sites,
            chrom_sizes = select_first([chrsz]),
            ligation_site = select_first([ligation_site]),
            quality = qualities[i],
            single_ended = align.bam_and_ligation_count[0][0].single_ended,
            runtime_environment = runtime_environment,
        }

        # If Juicer Tools doesn't support the assembly then need to pass chrom sizes
        if (!normalize_assembly_name.assembly_is_supported) {
            call create_hic as create_hic_with_chrom_sizes { input:
                pre = bam_to_pre.pre,
                pre_index = bam_to_pre.index,
                chrsz = select_first([chrsz]),
                restriction_sites = restriction_sites,
                quality = qualities[i],
                stats = calculate_stats.stats,
                stats_hists = calculate_stats.stats_hists,
                resolutions = create_hic_resolutions,
                assembly_name = assembly_name,
                num_cpus = create_hic_num_cpus,
                ram_gb = create_hic_ram_gb,
                juicer_tools_heap_size_gb = create_hic_juicer_tools_heap_size_gb,
                disk_size_gb = create_hic_disk_size_gb,
                runtime_environment = runtime_environment,
            }
        }

        if (normalize_assembly_name.assembly_is_supported) {
            call create_hic { input:
                pre = bam_to_pre.pre,
                pre_index = bam_to_pre.index,
                restriction_sites = restriction_sites,
                quality = qualities[i],
                stats = calculate_stats.stats,
                stats_hists = calculate_stats.stats_hists,
                resolutions = create_hic_resolutions,
                assembly_name = normalize_assembly_name.normalized_assembly_name,
                num_cpus = create_hic_num_cpus,
                ram_gb = create_hic_ram_gb,
                juicer_tools_heap_size_gb = create_hic_juicer_tools_heap_size_gb,
                disk_size_gb = create_hic_disk_size_gb,
                runtime_environment = runtime_environment,
            }
        }

        File unnormalized_hic_file = select_first([
            if (defined(create_hic.output_hic))
            then create_hic.output_hic
            else create_hic_with_chrom_sizes.output_hic
        ])

        call add_norm { input:
            hic = unnormalized_hic_file,
            normalization_methods = normalization_methods,
            quality = qualities[i],
            num_cpus = add_norm_num_cpus,
            ram_gb = add_norm_ram_gb,
            disk_size_gb = add_norm_disk_size_gb,
            runtime_environment = runtime_environment,
        }

        if (!no_call_tads) {
            call arrowhead { input:
                hic_file = add_norm.output_hic,
                quality = qualities[i],
                runtime_environment = runtime_environment,
            }
        }
        if (!no_call_loops) {
            if (!intact) {
                call hiccups { input:
                    hic_file = add_norm.output_hic,
                    quality = qualities[i],
                    runtime_environment = hiccups_runtime_environment,
                }
            }

            if (intact) {
                call hiccups_2 { input:
                    hic = add_norm.output_hic,
                    quality = qualities[i],
                    runtime_environment = hiccups_runtime_environment,
                }

                call localizer as localizer_intact { input:
                    hic = add_norm.output_hic,
                    loops = hiccups_2.merged_loops,
                    quality = qualities[i],
                    runtime_environment = runtime_environment,
                }
            }
        }

        if (defined(chrsz) && !no_eigenvectors) {
            call create_eigenvector { input:
                hic_file = add_norm.output_hic,
                chrom_sizes = select_first([chrsz]),
                output_filename_suffix = "_" + qualities[i],
                runtime_environment = runtime_environment,
            }

            call create_eigenvector as create_eigenvector_10kb { input:
                hic_file = add_norm.output_hic,
                chrom_sizes = select_first([chrsz]),
                resolution = 10000,
                output_filename_suffix = "_" + qualities[i],
                runtime_environment = runtime_environment,
            }
        }
    }

    if (!no_delta) {
        call delta { input:
            # Only run delta on MAPQ >= 30
            hic = if length(add_norm.output_hic) > 1 then add_norm.output_hic[1] else select_first([input_hic]),
            resolutions = delta_resolutions,
            models_path = delta_models_path,
            runtime_environment = delta_runtime_environment,
        }

        call localizer as localizer_delta { input:
            hic = if length(add_norm.output_hic) > 1 then add_norm.output_hic[1] else select_first([input_hic]),
            loops = delta.loops,
            runtime_environment = runtime_environment,
        }
    }

    if (defined(input_hic)) {
        if (!no_call_tads) {
            call arrowhead as arrowhead_input_hic { input:
                hic_file = select_first([input_hic]),
                runtime_environment = runtime_environment,
            }
        }
        if (!no_call_loops) {
            if (!intact) {
                call hiccups as hiccups_input_hic { input:
                    hic_file = select_first([input_hic]),
                    runtime_environment = hiccups_runtime_environment,
                }
            }

            if (intact) {
                call hiccups_2 as hiccups_2_input_hic { input:
                    hic = select_first([input_hic]),
                    runtime_environment = hiccups_runtime_environment,
                }
                call localizer as localizer_intact_input_hic { input:
                    hic = select_first([input_hic]),
                    loops = hiccups_2_input_hic.merged_loops,
                    runtime_environment = runtime_environment,
                }
            }
        }
    }

    if (!no_slice) {
        File hic_file = if length(add_norm.output_hic) > 1 then add_norm.output_hic[1] else select_first([input_hic])

        call slice as slice_25kb { input:
            hic_file = hic_file,
            resolution = 25000,
            runtime_environment = runtime_environment,
        }

        call slice as slice_50kb { input:
            hic_file = hic_file,
            resolution = 50000,
            runtime_environment = runtime_environment,
        }

        call slice as slice_100kb { input:
            hic_file = hic_file,
            resolution = 100000,
            runtime_environment = runtime_environment,
        }
    }
}

task get_ligation_site_regex {
    input {
        Array[String] restriction_enzymes
        RuntimeEnvironment runtime_environment
    }

    String output_path = "ligation_site_regex.txt"

    command <<<
        set -euo pipefail
        python3 "$(which get_ligation_site_regex.py)" \
            --enzymes ~{sep=" " restriction_enzymes} \
            --outfile ~{output_path}
    >>>

    output {
        String ligation_site_regex = read_string("~{output_path}")
        # Surface the original file for testing purposes
        File ligation_site_regex_file = "~{output_path}"
    }

    runtime {
        cpu : "1"
        memory: "500 MB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task normalize_assembly_name {
    input {
        String assembly_name
        String normalized_assembly_name_output_path = "normalized_assembly_name.txt"
        String assembly_is_supported_output_path = "is_supported.txt"
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        python3 "$(which normalize_assembly_name.py)" \
            ~{assembly_name} \
            ~{normalized_assembly_name_output_path} \
            ~{assembly_is_supported_output_path}
    >>>

    output {
        String normalized_assembly_name = read_string("~{normalized_assembly_name_output_path}")
        Boolean assembly_is_supported = read_boolean("~{assembly_is_supported_output_path}")
    }

    runtime {
        cpu : "1"
        memory: "500 MB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task align {
    input {
        FastqPair fastq_pair
        File idx_tar        # reference bwa index tar
        String ligation_site
        Int num_cpus = 32
        Int ram_gb = 64
        Int disk_size_gb = 1000
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        echo "Starting align"
        mkdir reference
        cd reference && tar -xvf ${idx_tar}
        index_folder=$(ls)
        reference_fasta=$(ls | head -1)
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ..

        usegzip=1
        name="result"
        ligation="${ligation_site}"
        name1=${fastq_pair.read_1}
        name2=${fastq_pair.read_2}
        ext=""
        singleend=~{if(defined(fastq_pair.read_2)) then "0" else "1"}
        #count ligations
        # Need to unset the -e option, when ligation site is XXXX grep will exit with
        # non-zero status
        set +e
        source /opt/scripts/common/countligations.sh
        set -e
        # Align reads
        echo "Running bwa command"
        bwa \
            mem \
            ~{if defined(fastq_pair.read_2) then "-SP5M" else "-5M"} \
            ~{if defined(fastq_pair.read_group) then "-R '" + fastq_pair.read_group + "'" else ""} \
            -t ~{num_cpus} \
            -K 320000000 \
            $reference_index_path \
            ${fastq_pair.read_1} \
            ~{default="" fastq_pair.read_2} | \
            samtools view -hbS - > aligned.bam
        mv result_norm.txt.res.txt ligation_count.txt
    }

    output {
        BamAndLigationCount bam_and_ligation_count = object {
            bam: "aligned.bam",
            ligation_count: "ligation_count.txt",
            single_ended: length(select_all([fastq_pair.read_2])) == 0,
        }
     }

    runtime {
        cpu : "~{num_cpus}"
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task chimeric_sam_specific {
    input {
        File bam
        File ligation_count
        File restriction_sites
        Boolean single_ended
        Int num_cpus = 8
        Int ram_gb = 16
        Int disk_size_gb = 1000
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        RESTRICTION_SITES_FILENAME=restriction_sites.txt
        gzip -dc ~{restriction_sites} > $RESTRICTION_SITES_FILENAME
        cp ~{ligation_count} result_norm.txt.res.txt
        samtools view -h -@ ~{num_cpus - 1} ~{bam} > result.sam
        awk \
            -v stem=result_norm \
            -v site_file=$RESTRICTION_SITES_FILENAME \
            ~{if(single_ended) then "-v singleend=1" else ""} \
            -f "$(which chimeric_sam.awk)" \
            result.sam | \
            samtools sort -t cb -n --threads ~{num_cpus} > chimeric_sam_specific.bam
    >>>

    output {
        File output_bam = "chimeric_sam_specific.bam"
        File stats = "result_norm.txt.res.txt"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory: "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task chimeric_sam_nonspecific {
    input {
        File bam
        File ligation_count
        Boolean single_ended
        Int num_cpus = 8
        Int ram_gb = 16
        Int disk_size_gb = 1000
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        cp ~{ligation_count} result_norm.txt.res.txt
        samtools view -h -@ ~{num_cpus - 1} ~{bam} > result.sam
        awk \
            -v stem=result_norm \
            ~{if(single_ended) then "-v singleend=1" else ""} \
            -f "$(which chimeric_sam.awk)" \
            result.sam > result.sam2
        ~{if(single_ended) then "samtools sort -t cb -n --threads " + num_cpus + " result.sam2 > chimeric_sam_nonspecific.bam && exit 0" else ""}
        awk \
            -v avgInsertFile=result_norm.txt.res.txt \
            -f "$(which adjust_insert_size.awk)" \
            result.sam2 | \
            samtools sort -t cb -n --threads ~{num_cpus} > chimeric_sam_nonspecific.bam
    >>>

    output {
        File output_bam = "chimeric_sam_nonspecific.bam"
        File stats = "result_norm.txt.res.txt"
    }

    runtime {
        cpu : "~{num_cpus}"
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task merge {
    input {
        Array[File] bams
        Int num_cpus = 8
        Int ram_gb = 16
        Int disk_size_gb = 6000
        String output_bam_filename = "merged"
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        samtools merge \
            -c \
            -t cb \
            -n \
            --threads ~{num_cpus - 1} \
            ~{output_bam_filename}.bam \
            ~{sep=' ' bams}
    >>>

    output {
        File bam = "~{output_bam_filename}.bam"
    }

    runtime {
        cpu : "~{num_cpus}"
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task dedup {
    input {
        File bam
        Int num_cpus = 8
        Int ram_gb = 32
        Int disk_size_gb = 5000
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        samtools view \
            -h \
            -@ ~{num_cpus - 1} \
            ~{bam} | \
            awk -f "$(which dups_sam.awk)" > merged_dedup.sam
        samtools view -b -@ ~{num_cpus - 1} merged_dedup.sam > merged_dedup.bam
    >>>

    output {
        File deduped_bam = "merged_dedup.bam"
    }

    runtime {
        cpu : "~{num_cpus}"
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task pre_to_pairs {
    input {
        File pre
        File chrom_sizes
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        PRE_FILENAME=pre.txt
        gzip -dc ~{pre} > $PRE_FILENAME
        perl "$(which juicer_shortform2pairs.pl)" $PRE_FILENAME ~{chrom_sizes} pairix
    }

    output {
        File out_file = "pairix.bsorted.pairs.gz"
    }

    runtime {
        cpu : "8"
        memory: "16 GB"
        disks: "local-disk 3000 HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task bam_to_pre {
    input {
        File bam
        Int quality
        Int num_cpus = 8
        String output_filename_suffix = ""
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        MERGED_NODUPS_FILENAME=merged_nodups_~{quality}~{output_filename_suffix}.txt
        MERGED_NODUPS_INDEX_FILENAME=merged_nodups_~{quality}~{output_filename_suffix}_index.txt
        samtools view \
            -h \
            -F 1024 \
            -O sam \
            ~{bam} \
            -@ ~{num_cpus - 1} | \
            awk -v mapq=~{quality} -f "$(which sam_to_pre.awk)" > $MERGED_NODUPS_FILENAME
        $(which index_by_chr.awk) $MERGED_NODUPS_FILENAME 500000 > $MERGED_NODUPS_INDEX_FILENAME
        gzip -n $MERGED_NODUPS_FILENAME
        gzip -n $MERGED_NODUPS_INDEX_FILENAME
    >>>

    output {
        File pre = "merged_nodups_~{quality}~{output_filename_suffix}.txt.gz"
        File index = "merged_nodups_~{quality}~{output_filename_suffix}_index.txt.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 3000 HDD"
        memory : "64 GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task calculate_stats {
    input {
        Array[File] alignment_stats
        File pre
        File chrom_sizes
        File bam
        File? restriction_sites
        String ligation_site
        String output_filename_suffix = ""
        Boolean single_ended = false
        Int quality = 0
        RuntimeEnvironment runtime_environment
    }

    command <<<
        PRE_FILE=pre.txt
        RESTRICTION_SITES_FILENAME=restriction_sites.txt
        STATS_FILENAME=stats_~{quality}~{output_filename_suffix}.txt
        gzip -dc ~{pre} > $PRE_FILE
        ~{if defined(restriction_sites) then "gzip -dc " + restriction_sites + " > $RESTRICTION_SITES_FILENAME" else ""}
        if [ ~{if(single_ended) then "1" else "0"} -eq 1 ]
        then
            RET=$(samtools view -f 1024 -F 256 ~{bam} | awk '{if ($0~/rt:A:7/){singdup++}else{dup++}}END{print dup,singdup}')
            DUPS=$(echo $RET | awk '{print $1}')
            SINGDUPS=$(echo $RET | awk '{print $2}')
            awk \
                -f "$(which stats_sub.awk)" \
                -v dups=$DUPS \
                -v singdups=$SINGDUPS \
                -v ligation="~{ligation_site}" \
                -v singleend=1 \
                ~{sep=" " alignment_stats} >> $STATS_FILENAME
        else
            DUPS=$(samtools view -c -f 1089 -F 256 ~{bam})
            awk \
                -f "$(which stats_sub.awk)" \
                -v dups=$DUPS \
                -v ligation=$ligation \
                ~{sep=" " alignment_stats} >> $STATS_FILENAME
        fi
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx16g \
            -jar /opt/scripts/common/juicer_tools.jar \
            statistics \
            ~{if defined(restriction_sites) then "$RESTRICTION_SITES_FILENAME" else "none"} \
            $STATS_FILENAME \
            ~{pre} \
            ~{chrom_sizes}
        python3 \
            "$(which jsonify_stats.py)" \
            $STATS_FILENAME \
            stats_~{quality}~{output_filename_suffix}.json
    >>>

    output {
        File stats = "stats_~{quality}~{output_filename_suffix}.txt"
        File stats_json = "stats_~{quality}~{output_filename_suffix}.json"
        File stats_hists = "stats_~{quality}~{output_filename_suffix}_hists.m"
    }

    runtime {
        cpu : "1"
        disks: "local-disk 4000 HDD"
        memory : "16 GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task create_hic {
    input {
        File pre
        File pre_index
        File stats
        File stats_hists
        Int quality
        Array[Int] resolutions
        String? assembly_name
        File? chrsz
        File? restriction_sites
        Int num_cpus = 16
        # Should always set juicer_tools_heap_size_gb < ram_gb
        Int ram_gb = 384
        Int juicer_tools_heap_size_gb = 370
        Int disk_size_gb = 2000
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        PRE_FILE=pre.txt
        PRE_INDEX_FILE=pre_index.txt
        RESTRICTION_SITES_FILENAME=restriction_sites.txt
        gzip -dc ~{pre} > $PRE_FILE
        gzip -dc ~{pre_index} > $PRE_INDEX_FILE
        ~{if defined(restriction_sites) then "gzip -dc " + restriction_sites + " > $RESTRICTION_SITES_FILENAME" else ""}
        # If the assembly name is empty, then we write chrsz path into file as usual, otherwise, use the assembly name instead of the path
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx~{juicer_tools_heap_size_gb}g \
            -jar /opt/scripts/common/juicer_tools.jar \
            pre \
            -n \
            ~{if defined(restriction_sites) then "-f $RESTRICTION_SITES_FILENAME" else ""} \
            -s ~{stats} \
            -g ~{stats_hists} \
            ~{if defined(assembly_name) then "-y " + assembly_name else ""} \
            -r ~{sep="," resolutions} \
            -i $PRE_INDEX_FILE \
            --block-capacity 1000000 \
            --threads ~{num_cpus} \
            $PRE_FILE \
            inter_~{quality}_unnormalized.hic \
            ~{if defined(chrsz) then chrsz else assembly_name}
    >>>

    output {
        File output_hic = "inter_~{quality}_unnormalized.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        memory : "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}


task add_norm {
    input {
        File hic
        Array[String] normalization_methods = []
        Int? mthreads #only used with addnorm2
        Int? save_ram #only used with addnorm2
        Int quality
        Int num_cpus = 24
        Int disk_size_gb = 256
        Int ram_gb = 72
        String juicer_tools_jar = "/opt/scripts/common/juicer_tools.jar"
        String normalization_command = "addNorm"
        RuntimeEnvironment runtime_environment
    }

    Int java_heap_gb = ram_gb - 12

    command {
        set -euo pipefail
        cp ~{hic} inter_~{quality}.hic
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx~{if java_heap_gb > 0 then java_heap_gb else 12}g \
            -jar ~{juicer_tools_jar} \
            ~{normalization_command} \
            ~{if length(normalization_methods) > 0 then "-k" else ""} ~{sep="," normalization_methods} \
            ~{"--mthreads" + mthreads} \
            ~{"--save-ram" + save_ram} \
            --threads ~{num_cpus} \
            inter_~{quality}.hic
    }

    output {
        File output_hic = "inter_~{quality}.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory : "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}


task arrowhead {
    input {
        File hic_file
        Int disk_size_gb = 100
        Int num_cpus = 1
        Int ram_gb = 32
        Int quality = 0
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx16g \
            -jar /opt/scripts/common/juicer_tools.jar \
            arrowhead \
            ~{hic_file} \
            contact_domains
        gzip -n contact_domains/*
        STEM=$(basename contact_domains/*.bedpe.gz .bedpe.gz)
        mv contact_domains/*.bedpe.gz "${STEM}_~{quality}.bedpe.gz"
    >>>

    output {
        File out_file = glob('*_~{quality}.bedpe.gz')[0]
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory : "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task hiccups {
    input {
        File hic_file
        Int quality = 0
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -jar /opt/scripts/common/juicer_tools.jar \
            hiccups \
            ${hic_file} \
            loops
        gzip -n loops/merged_loops.bedpe
        mv loops/merged_loops.bedpe.gz merged_loops_~{quality}.bedpe.gz
    }

    output {
        File merged_loops = "merged_loops_~{quality}.bedpe.gz"
    }

    runtime {
        cpu : "1"
        bootDiskSizeGb: "20"
        disks: "local-disk 100 HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        memory: "8 GB"
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

task hiccups_2 {
    input {
        File hic
        Int quality = 0
        Int disk_size_gb = 100
        Int num_cpus = 2
        Int num_gpus = 1
        Int ram_gb = 64
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx60G \
            -Xms60G \
            -jar /opt/feature_tools.jar \
            hiccups2 \
            -k SCALE \
            --threads ~{num_cpus} \
            ~{hic} \
            loops
        gzip -n loops/merged_loops.bedpe
        mv loops/merged_loops.bedpe.gz merged_loops_~{quality}.bedpe.gz
    }

    output {
        File merged_loops = "merged_loops_~{quality}.bedpe.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        bootDiskSizeGb: "20"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
        gpuType: "nvidia-tesla-p100"
        gpuCount: "~{num_gpus}"
        memory: "~{ram_gb} GB"
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

task localizer {
    input {
        File hic
        File loops
        Int localizer_resolution = 100
        Int localizer_window = 1
        Int quality = 0
        Int disk_size_gb = 100
        Int num_cpus = 24
        Int ram_gb = 64
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        export LOOPS_FILE=loops.bedpe
        gzip -dc ~{loops} > $LOOPS_FILE
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx60G \
            -Xms60G \
            -jar /opt/feature_tools.jar \
            localizer \
            -k SCALE \
            -r ~{localizer_resolution} \
            -w ~{localizer_window} \
            --threads ~{num_cpus} \
            ~{hic} \
            $LOOPS_FILE \
            localized
        gzip -n localized/localizedList_primary_~{localizer_resolution}.bedpe
        mv localized/localizedList_primary_~{localizer_resolution}.bedpe.gz localized_loops_~{quality}.bedpe.gz
    }

    output {
        File localized_loops = "localized_loops_~{quality}.bedpe.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory: "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task delta {
    input {
        File hic
        Array[Int] resolutions
        Float threshold = 0.85
        Int disk_size_gb = 100
        Int gpu_count = 1
        Int ram_gb = 32
        String normalization = "SCALE"
        String stem = "predicted"
        String models_path = "beta-models"
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        python \
            "$(which Deploy.py)" \
            ~{hic} \
            /opt/deploy-delta/~{models_path} \
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
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
        gpuType: "nvidia-tesla-p100"
        gpuCount: "~{gpu_count}"
        memory: "~{ram_gb} GB"
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

task create_eigenvector {
    input {
        File hic_file
        File chrom_sizes
        Int disk_size_gb = 100
        Int ram_gb = 8
        Int num_cpus = 16
        Int resolution = 5000
        String output_filename_suffix = ""
        String normalization = "SCALE"
        RuntimeEnvironment runtime_environment
    }

    command {
        newGW_Intra_Flip \
            -n ~{normalization} \
            -T ~{num_cpus} \
            ~{hic_file} \
            eigenvector_~{resolution}~{output_filename_suffix}.wig \
            ~{resolution}
        wigToBigWig \
            eigenvector_~{resolution}~{output_filename_suffix}.wig \
            ~{chrom_sizes} \
            eigenvector_~{resolution}~{output_filename_suffix}.bw
    }

    output {
        File eigenvector_wig = "eigenvector_~{resolution}~{output_filename_suffix}.wig"
        File eigenvector_bigwig = "eigenvector_~{resolution}~{output_filename_suffix}.bw"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory : "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task slice {
    input {
        File hic_file
        Int disk_size_gb = 100
        Int ram_gb = 24
        Int num_cpus = 1
        Int resolution = 25000
        Int minimum_num_clusters = 2
        Int maximum_num_clusters = 13
        Int num_kmeans_runs = 4
        RuntimeEnvironment runtime_environment
    }

    command {
        set -euo pipefail
        java \
            -Xmx20G \
            -jar /opt/MixerTools.jar \
            slice \
            --encode-mode \
            -r ~{resolution} \
            ~{hic_file} \
            ~{minimum_num_clusters},~{maximum_num_clusters},~{num_kmeans_runs} \
            slice_results \
            cell_type
        gzip -n slice_results/*.bed
        mv slice_results/slice_subcompartment_clusters.bed.gz slice_subcompartment_clusters_~{resolution}.bed.gz
    }

    output {
        File subcompartments = "slice_subcompartment_clusters_~{resolution}.bed.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} SSD"
        memory : "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task create_accessibility_track {
    input {
        File pre
        File chrom_sizes
        Int ram_gb = 64
        Int disk_size_gb = 1000
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        PRE_FILE=pre.txt
        gzip -dc ~{pre} > $PRE_FILE
        awk '{print $1}' ~{chrom_sizes} | while read chrom; do awk -v chr=${chrom} 'BEGIN{OFS="\t"}$2==chr{c[$3]++}$6==chr{c[$7]++}END{for (i in c) {print chr, i-1, i, c[i]}}' $PRE_FILE | sort -k2,2n >> merged30.bedgraph; done;
        sort -k1,1 -k2,2n -S6G merged30.bedgraph > merged30.sorted.bedgraph
        bedGraphToBigWig merged30.sorted.bedgraph ~{chrom_sizes} inter_30.bw
    >>>

    output {
        File bigwig = "inter_30.bw"
    }

    runtime {
        cpu : "1"
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task exit_early {
    input {
        String message
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        echo ~{message}
        exit 1
    >>>

    runtime {
        cpu : "1"
        memory: "500 MB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
