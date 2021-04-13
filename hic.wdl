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
        version: "0.1.0"
        caper_docker: "gcr.io/hic-pipeline/new-jar:latest"
        caper_singularity: "docker://encodedcc/hic-pipeline:0.1.0"
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
        String? assembly_name
    }

    parameter_meta {
        fastq: "Twice nested array of input `FastqPair`s, takes form of [lib_id][fastq_id]"
        restriction_enzymes: "An array of names containing the restriction enzyme(s) used to generate the Hi-C libraries"
        restriction_sites: "A text file containing cut sites for the given restriction enzyme. You should generate this file using this script: https://github.com/aidenlab/juicer/blob/encode/misc/generate_site_positions.py"
        ligation_site_regex: "A custom regex to use for counting ligation site, if specified then restriction_sites file must be manually specified. Can be just a single site, e.g. ATGC, or several sites wrapped in parentheses and separated by pipes, e.g. `(ATGC|CTAG)`"
        chrsz: "A chromosome sizes file for the desired assembly, this is a tab-separated text file whose rows take the form [chromosome] [size]"
        reference_index: "A pregenerated BWA index for the desired assembly"
        normalization_methods: "An array of normalization methods to use for .hic file generation as per Juicer Tools `pre`, if not specified then will use `pre` defaults of VC, VC_SQRT, KR, and SCALE. Valid methods are VC, VC_SQRT, KR, SCALE, GW_KR, GW_SCALE, GW_VC, INTER_KR, INTER_SCALE, and INTER_VC."
        input_pairs: "A text file containing the paired fragments to use to generate the .hic contact maps, a detailed format description can be found here: https://github.com/aidenlab/juicer/wiki/Pre#long-format"
        input_pairs_index: "Index of input_pairs as generated with index_by_chr.awk in task bam_to_pre"
        input_hic: "An input .hic file for which to call loops and domains"
        no_bam2pairs: "If set to `true`, avoid generating .pairs files, defaults to false"
        no_call_loops: "If set to `true`, avoid calling loops with hiccups, defaults to false"
        no_call_tads: "If set to `true`, avoid calling domains with arrowhead, defaults to false"
        align_num_cpus: "Number of threads to use for bwa alignment"
        assembly_name: "Name of assembly to insert into hic file header, recommended to specify for reproducbility otherwise hic file will be nondeterministic"
    }

    # Default MAPQ thresholds for generating .hic contact maps
    Array[Int] DEFAULT_HIC_QUALITIES = [1, 30]
    Boolean is_nonspecific = length(restriction_enzymes) > 0 && restriction_enzymes[0] == "none"

    if (!defined(ligation_site_regex) && !defined(input_hic)) {
        call get_ligation_site_regex { input:
            restriction_enzymes = restriction_enzymes
        }
        String ligation_site = select_first([ligation_site_regex, get_ligation_site_regex.ligation_site_regex])

        if (!is_nonspecific && !defined(restriction_sites)) {
            call exit_early { input:
                message = "Must provide restriction sites file if enzyme is not `none`"
            }
        }
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
            duplicate_counts = dedup.duplicate_count,
            pre = bam_to_pre.pre,
            restriction_sites = restriction_sites,
            ligation_site = select_first([ligation_site]),
            quality = qualities[i],
        }

        call create_hic { input:
            pre = select_first([input_pairs, bam_to_pre.pre]),
            pre_index = select_first([input_pairs_index, bam_to_pre.index]),
            chrsz = select_first([chrsz]),
            quality = qualities[i],
            stats = calculate_stats.stats,
            stats_hists = calculate_stats.stats_hists,
            assembly_name = assembly_name,
            normalization_methods = normalization_methods,
        }
    }

    if ((defined(input_hic) || defined(create_hic.output_hic))) {
        File hic_file = if defined(input_hic) then select_first([input_hic]) else create_hic.output_hic[1]
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
        cp ~{ligation_count} result_norm.txt.res.txt
        samtools view -h -@ ~{num_cpus - 1} ~{bam} > result.sam
        awk \
            -v stem=result_norm \
            -v site_file=~{restriction_sites} \
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
        disks: "local-disk 1000 HDD"
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
            awk -f "$(which dups_sam.awk)" -v fname=duplicate_count.txt > merged_dedup.sam
        samtools view -b -@ ~{num_cpus - 1} merged_dedup.sam > merged_dedup.bam
    >>>

    output {
        File deduped_bam = "merged_dedup.bam"
        File duplicate_count = "duplicate_count.txt"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 2000 SSD"
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
        disks: "local-disk 1000 HDD"
    }
}

task bam_to_pre {
    input {
        File bam
        Int quality
        Int num_cpus = 8
    }

    command <<<
        set -euo pipefail
        MERGED_NODUPS_FILENAME=merged_nodups_~{quality}.txt
        MERGED_NODUPS_INDEX_FILENAME=merged_nodups_~{quality}_index.txt
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
        File pre = "merged_nodups_~{quality}.txt.gz"
        File index = "merged_nodups_~{quality}_index.txt.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
    }
}

task calculate_stats {
    input {
        Array[File] alignment_stats
        Array[File] duplicate_counts
        File pre
        File? restriction_sites
        String ligation_site
        Int quality
    }

    command <<<
        PRE_FILE=pre.txt
        STATS_FILENAME=stats_~{quality}.txt
        gzip -dc ~{pre} > $PRE_FILE
        awk -f "$(which stats_sub.awk)" ~{sep=" " alignment_stats} >> $STATS_FILENAME
        awk \
            -f "$(which count_unique_reads.awk)" \
            $STATS_FILENAME \
            ~{sep=" " duplicate_counts} >> $STATS_FILENAME
        statistics.pl \
            -s ~{default="none" restriction_sites} \
            -l ~{ligation_site} \
            -o $STATS_FILENAME \
            $PRE_FILE
        python3 "$(which jsonify_stats.py)" --alignment-stats $STATS_FILENAME
    >>>

    output {
        File stats = "stats_~{quality}.txt"
        File stats_json = "stats_~{quality}.json"
        File stats_hists = "stats_~{quality}_hists.m"
    }

    runtime {
        cpu : "1"
        disks: "local-disk 1000 HDD"
        memory : "8 GB"
    }
}

task create_hic {
    input {
        File chrsz
        File pre
        File pre_index
        File stats
        File stats_hists
        Array[String] normalization_methods = []
        Int quality
        String? assembly_name
        Int num_cpus = 16
    }

    command <<<
        set -euo pipefail
        PRE_FILE=pre.txt
        PRE_INDEX_FILE=pre_index.txt
        gzip -dc ~{pre} > $PRE_FILE
        gzip -dc ~{pre_index} > $PRE_INDEX_FILE
        # If the assembly name is empty, then we write chrsz path into file as usual, otherwise, use the assembly name instead of the path
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx240g \
            -jar /opt/scripts/common/juicer_tools.jar \
            pre \
            -n \
            -s ~{stats} \
            -g ~{stats_hists} \
            -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 \
            ~{if defined(assembly_name) then "-y " + assembly_name else ""} \
            ~{if length(normalization_methods) > 0 then "-k" else ""} ~{sep="," normalization_methods} \
            -i $PRE_INDEX_FILE \
            --threads ~{num_cpus} \
            $PRE_FILE \
            inter_~{quality}.hic \
            ~{chrsz}
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx240g \
            -jar /opt/scripts/common/juicer_tools.jar \
            addNorm \
            --threads ~{num_cpus} \
            inter_~{quality}.hic
    >>>

    output {
        File output_hic = "inter_~{quality}.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 2000 SSD"
        memory : "256 GB"
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
        memory : "16 GB"
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
    }

    runtime {
        cpu : "1"
        bootDiskSizeGb: "20"
        disks: "local-disk 100 SSD"
        docker: "gcr.io/hic-pipeline/new-jar-hiccups:latest"
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
