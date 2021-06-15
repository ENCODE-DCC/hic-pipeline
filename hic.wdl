version 1.0

struct FastqPair {
    File read_1
    File read_2
    String? read_group
}

struct BamAndLigationCount {
    File bam
    File ligation_count
}

workflow hic {
    meta {
        version: "0.6.0"
        caper_docker: "encodedcc/hic-pipeline:0.6.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:0.6.0"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json"
    }

    input {
        # Main entrypoint, need to specify all five of these values except read_groups when running from fastqs
        Array[Array[FastqPair]] fastq = []
        Array[String] restriction_enzymes = []
        String? ligation_site_regex
        File? restriction_sites
        File? chrsz
        File? reference_index

        # Entrypoint right before hic generation
        File? input_pairs
        File? input_pairs_index

        # Entrypoint for loop and TAD calls
        File? input_hic

        Array[String] normalization_methods = []
        Boolean no_bam2pairs = false
        Boolean no_call_loops = false
        Boolean no_call_tads = false
        Int align_num_cpus = 32
        Int? create_hic_num_cpus
        String assembly_name = "undefined"
    }

    # Default MAPQ thresholds for generating .hic contact maps
    Array[Int] DEFAULT_HIC_QUALITIES = [1, 30]
    Boolean is_nonspecific = length(restriction_enzymes) > 0 && restriction_enzymes[0] == "none"

    if (!defined(input_hic)) {
        if (!defined(ligation_site_regex)) {
            call get_ligation_site_regex { input:
                restriction_enzymes = restriction_enzymes
            }
        }

        String ligation_site = select_first([ligation_site_regex, get_ligation_site_regex.ligation_site_regex])

        if (!is_nonspecific && !defined(restriction_sites)) {
            call exit_early { input:
                message = "Must provide restriction sites file if enzyme is not `none`"
            }
        }
    }

    call normalize_assembly_name { input:
        assembly_name = assembly_name
    }

    scatter(i in range(length(fastq))) {
        Array[FastqPair] replicate = fastq[i]
        scatter(fastq_pair in replicate) {
            call align { input:
                fastq_pair = fastq_pair,
                idx_tar = select_first([reference_index]),
                ligation_site = select_first([ligation_site]),
                num_cpus = align_num_cpus,
            }
        }

        if (is_nonspecific) {
            scatter(bam_and_ligation_count in align.bam_and_ligation_count) {
                call chimeric_sam_nonspecific { input:
                    bam = bam_and_ligation_count.bam,
                    ligation_count = bam_and_ligation_count.ligation_count,
                }
            }
        }

        if (!is_nonspecific) {
            scatter(bam_and_ligation_count in align.bam_and_ligation_count) {
                call chimeric_sam_specific { input:
                    bam = bam_and_ligation_count.bam,
                    ligation_count = bam_and_ligation_count.ligation_count,
                    restriction_sites = select_first([restriction_sites]),
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
            output_bam_filename = "merged_" + i
        }

        call dedup { input:
            bam = merge.bam
        }

        call bam_to_pre as bam_to_pre_for_stats { input:
            bam = dedup.deduped_bam,
            quality = 1,
            output_filename_suffix = "_lib" + i
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
        }
    }

    if (!defined(input_hic)) {
        call merge as merge_replicates { input:
            bams = dedup.deduped_bam,
        }
        # convert alignable bam to pairs to be consistent with 4DN
        if ( !no_bam2pairs && defined(chrsz)) {
            call bam2pairs { input:
                bam_file = merge_replicates.bam,
                chrsz_  = select_first([chrsz])
            }
        }
    }

    Array[Int] qualities = if !defined(input_hic) then DEFAULT_HIC_QUALITIES else []
    scatter(i in range(length(qualities))) {
        call bam_to_pre { input:
            bam = select_first([merge_replicates.bam]),
            quality = qualities[i],
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
        }

        # If Juicer Tools doesn't support the assembly then need to pass chrom sizes
        if (!normalize_assembly_name.assembly_is_supported) {
            call create_hic as create_hic_with_chrom_sizes { input:
                pre = select_first([input_pairs, bam_to_pre.pre]),
                pre_index = select_first([input_pairs_index, bam_to_pre.index]),
                chrsz = select_first([chrsz]),
                restriction_sites = restriction_sites,
                quality = qualities[i],
                stats = calculate_stats.stats,
                stats_hists = calculate_stats.stats_hists,
                assembly_name = assembly_name,
                normalization_methods = normalization_methods,
                num_cpus = create_hic_num_cpus,
            }
        }

        if (normalize_assembly_name.assembly_is_supported) {
            call create_hic { input:
                pre = select_first([input_pairs, bam_to_pre.pre]),
                pre_index = select_first([input_pairs_index, bam_to_pre.index]),
                restriction_sites = restriction_sites,
                quality = qualities[i],
                stats = calculate_stats.stats,
                stats_hists = calculate_stats.stats_hists,
                assembly_name = normalize_assembly_name.normalized_assembly_name,
                normalization_methods = normalization_methods,
                num_cpus = create_hic_num_cpus,
            }
        }
    }

    File hic_file = select_first(
        [input_hic, create_hic.output_hic[1], create_hic_with_chrom_sizes.output_hic[1]]
    )
    if (!no_call_tads) {
        call arrowhead { input:
            hic_file = hic_file
        }
    }
    if (!no_call_loops) {
        call hiccups { input:
            hic_file = hic_file
        }
    }
}

task get_ligation_site_regex {
    input {
        Array[String] restriction_enzymes
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
    }
}

task normalize_assembly_name {
    input {
        String assembly_name
        String normalized_assembly_name_output_path = "normalized_assembly_name.txt"
        String assembly_is_supported_output_path = "is_supported.txt"
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
    }
}

task align {
    input {
        FastqPair fastq_pair
        File idx_tar        # reference bwa index tar
        String ligation_site
        Int num_cpus = 32
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
            -SP5M \
            ~{if defined(fastq_pair.read_group) then "-R '" + fastq_pair.read_group + "'" else ""} \
            -t ~{num_cpus} \
            -K 320000000 \
            $reference_index_path \
            ${fastq_pair.read_1} \
            ${fastq_pair.read_2} | \
            samtools view -hbS - > aligned.bam
        mv result_norm.txt.res.txt ligation_count.txt
    }

    output {
        BamAndLigationCount bam_and_ligation_count = object {
            bam: "aligned.bam",
            ligation_count: "ligation_count.txt",
        }
     }

    runtime {
        cpu : "~{num_cpus}"
        memory: "64 GB"
        disks: "local-disk 1000 HDD"
    }
}

task chimeric_sam_specific {
    input {
        File bam
        File ligation_count
        File restriction_sites
        Int num_cpus = 8
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
        disks: "local-disk 1000 HDD"
        memory: "16 GB"
    }
}

task chimeric_sam_nonspecific {
    input {
        File bam
        File ligation_count
        Int num_cpus = 8
    }

    command <<<
        set -euo pipefail
        cp ~{ligation_count} result_norm.txt.res.txt
        samtools view -h -@ ~{num_cpus - 1} ~{bam} > result.sam
        awk \
            -v stem=result_norm \
            -f "$(which chimeric_sam.awk)" \
            result.sam > result.sam2
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
        disks: "local-disk 1000 HDD"
        memory: "16 GB"
    }
}

task merge {
    input {
        Array[File] bams
        Int num_cpus = 8
        String output_bam_filename = "merged"
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
        memory: "16 GB"
        disks: "local-disk 6000 HDD"
    }
}

task dedup {
    input {
        File bam
        Int num_cpus = 8
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
        disks: "local-disk 5000 SSD"
        memory: "32 GB"
    }
}

task bam2pairs {
    input {
        File bam_file
        File chrsz_
    }

    command {
        set -euo pipefail
        bam2pairs -c ${chrsz_} ${bam_file} pairix
    }

    output {
        File out_file = "pairix.bsorted.pairs.gz"
        File pairs_index = "pairix.bsorted.pairs.gz.px2"
    }

    runtime {
        cpu : "8"
        memory: "16 GB"
        disks: "local-disk 3000 HDD"
    }
}

task bam_to_pre {
    input {
        File bam
        Int quality
        Int num_cpus = 8
        String output_filename_suffix = ""
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
        Int quality = 0
    }

    command <<<
        PRE_FILE=pre.txt
        RESTRICTION_SITES_FILENAME=restriction_sites.txt
        STATS_FILENAME=stats_~{quality}~{output_filename_suffix}.txt
        gzip -dc ~{pre} > $PRE_FILE
        ~{if defined(restriction_sites) then "gzip -dc " + restriction_sites + " > $RESTRICTION_SITES_FILENAME" else ""}
        duplicate_count=$(samtools view -c -f 1089 -F 256 ~{bam})
        awk \
            -f "$(which stats_sub.awk)" \
            -v ligation="~{ligation_site}" \
            -v dups="$duplicate_count" \
            ~{sep=" " alignment_stats} >> $STATS_FILENAME
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx16g \
            -jar /opt/scripts/common/juicer_tools.jar \
            statistics \
            --ligation "~{ligation_site}" \
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
    }
}

task create_hic {
    input {
        File pre
        File pre_index
        File stats
        File stats_hists
        Array[String] normalization_methods = []
        Int quality
        String? assembly_name
        File? chrsz
        File? restriction_sites
        Int? num_cpus = 8
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
            -Xmx590g \
            -jar /opt/scripts/common/juicer_tools.jar \
            pre \
            -n \
            ~{if defined(restriction_sites) then "-f $RESTRICTION_SITES_FILENAME" else ""} \
            -s ~{stats} \
            -g ~{stats_hists} \
            ~{if defined(assembly_name) then "-y " + assembly_name else ""} \
            -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 \
            -i $PRE_INDEX_FILE \
            --block-capacity 1000000 \
            --threads ~{num_cpus} \
            $PRE_FILE \
            inter_~{quality}.hic \
            ~{if defined(chrsz) then chrsz else assembly_name}
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx590g \
            -jar /opt/scripts/common/juicer_tools.jar \
            addNorm \
            ~{if length(normalization_methods) > 0 then "-k" else ""} ~{sep="," normalization_methods} \
            --threads ~{num_cpus} \
            inter_~{quality}.hic
    >>>

    output {
        File output_hic = "inter_~{quality}.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 2000 SSD"
        memory : "600 GB"
    }
}


task arrowhead {
    input {
        File hic_file
    }

    command {
        set -euo pipefail
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx16g \
            -jar /opt/scripts/common/juicer_tools.jar \
            arrowhead \
            ${hic_file} \
            contact_domains
        gzip -n contact_domains/*
    }

    output {
        File out_file = glob('contact_domains/*.bedpe.gz')[0]
    }

    runtime {
        cpu : "1"
        disks: "local-disk 100 SSD"
        memory : "32 GB"
    }
}

task hiccups {
    input {
        File hic_file
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
        gzip -n loops/*.bedpe
    }

    output {
        File out_file = glob("loops/*.bedpe.gz")[0]
        File merged_loops = "loops/merged_loops.bedpe.gz"
    }

    runtime {
        cpu : "1"
        bootDiskSizeGb: "20"
        disks: "local-disk 100 SSD"
        docker: "encodedcc/hic-pipeline:0.6.0_hiccups"
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

task exit_early {
    input {
        String message
    }

    command <<<
        set -euo pipefail
        echo ~{message}
        exit 1
    >>>

    runtime {
        cpu : "1"
        memory: "500 MB"
    }
}
