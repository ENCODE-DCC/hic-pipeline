version 1.0

#CAPER docker encodedcc/hic-pipeline:0.1.0
#CAPER singularity docker://encodedcc/hic-pipeline:0.1.0
#CROO out_def https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json

workflow hic {
    input {
        # Main entrypoint, need to specify all five of these values when running from fastqs
        Array[Array[Array[File]]] fastq = []
        Array[Array[String]] read_groups = []
        Array[String] restriction_enzymes
        File? restriction_sites
        File? chrsz
        File? reference_index

        # Entrypoint from aligned bam
        Array[Array[File]]? bams
        Array[Array[File]]? ligation_counts

        # Entrypoint right before hic generation
        File? input_pairs

        # Entrypoint for loop and TAD calls
        File? input_hic

        # Inputs for entrypoint after library processing
        Array[File]? input_dedup_pairs
        Array[File]? alignment_stats
        Array[File]? library_stats

        # Input to build restriction site locations
        File? reference_fasta
        Boolean restriction_site_locations_only = false

        Boolean no_bam2pairs = false
        Boolean no_call_loops = false
        Boolean no_call_tads = false
        Boolean include_mapq0_reads = false
        Int cpu = 32
        String? assembly_name
    }

    parameter_meta {
        fastq: "Twice nested array of input fastqs, takes form of [lib_id][fastq_id][read_end_id]"
        read_groups: "Optional strings to be inserted into the BAM as the read group (@RG), passed via `samtools addreplacerg` `-r` option. One per SE read/read pair with nested array structure mirroring the `fastq` input"
        restriction_enzyme: "An array of names containing the restriction enzyme(s) used to generate the Hi-C libraries"
        restriction_sites: "A text file containing cut sites for the given restriction enzyme. You should generate this file using this script: https://github.com/aidenlab/juicer/blob/encode/misc/generate_site_positions.py"
        chrsz: "A chromosome sizes file for the desired assembly, this is a tab-separated text file whose rows take the form [chromosome] [size]"
        reference_index: "A pregenerated BWA index for the desired assembly"
        bams: "Aligned, unfiltered bams, organized by [biorep[techrep]]. If specified, the `ligation_counts` array must also be specified"
        ligation_counts: "Text files containing ligation counts for the fastq pair, organized by [biorep[techrep]]. Has no meaning if the `bams` array is not also be specified. These should be calculated from fastqs using the Juicer countligations script: https://github.com/aidenlab/juicer/blob/encode/CPU/common/countligations.sh"
        input_pairs: "A text file containing the paired fragments to use to generate the .hic contact maps, a detailed format description can be found here: https://github.com/aidenlab/juicer/wiki/Pre#long-format"
        input_hic: "An input .hic file for which to call loops and domains"
        input_dedup_pairs: "An optional array consisting of text files of paired fragments, one per library, same format as input_pairs, used for merging libraries"
        alignment_stats: "An optional array consisting of text files of alignment stats, one per library, only has meaning when used in combination with `input_dedup_pairs`. Use is recommended but not required when merging libraries in order to calculate quality metrics on the merged libraries."
        library_stats: "An optional array consisting of text files of library stats, one per library, only has meaning when used in combination with `input_dedup_pairs`. Use is recommended but not required when merging libraries in order to calculate quality metrics on the merged libraries."
        reference_fasta: "FASTA file for the genome of interest to be used for generating restriction site locations. For the output locations file to have a descriptive filename it is also recommended to specify the `assembly_name`. Has no use if a pregenerated restriction site locations file is provided."
        restriction_site_locations_only: "If `true`, then will only generate the restriction site locations file."
        no_bam2pairs: "If set to `true`, avoid generating .pairs files, defaults to false"
        no_call_loops: "If set to `true`, avoid calling loops with hiccups, defaults to false"
        no_call_tads: "If set to `true`, avoid calling domains with arrowhead, defaults to false"
        include_mapq0_reads: "If set to `true`, chimeric reads (3 alignments) with one MAPQ 0 read will be classified as normal paired reads with the MAPQ 0 read discarded. If `false`, such reads will be classified as low mapq collisions, defaults to false"
        cpu: "Number of threads to use for bwa alignment"
        assembly_name: "Name of assembly to insert into hic file header, recommended to specify for reproducbility otherwise hic file will be nondeterministic"
    }

    # Default MAPQ thresholds for generating .hic contact maps
    Array[String] DEFAULT_HIC_QUALITIES = ["1", "30"]

    if (length(restriction_enzymes) > 1 && !defined(restriction_sites)) {
        call exit_early { input:
            message = "To use multiple restriction enzymes you must generate the restriction sites file manually"
        }
    }

    call get_ligation_site_regex { input:
        restriction_enzymes = restriction_enzymes
    }

    String ligation_site = get_ligation_site_regex.ligation_site_regex

    # make_restriction_site_locations currently supports only supports one enzyme
    if (defined(reference_fasta) && !defined(restriction_sites) && length(restriction_enzymes) == 1) {
        call make_restriction_site_locations { input:
            reference_fasta = select_first([reference_fasta]),
            assembly_name = select_first([assembly_name, "unknown_assembly"]),
            restriction_enzyme = restriction_enzymes[0],
        }
    }

    Boolean has_restriction_sites = defined(restriction_sites) || defined(make_restriction_site_locations.restriction_site_locations)

    # scatter over libraries
    if (has_restriction_sites && defined(chrsz) && !restriction_site_locations_only) {
        File chrsz_ = select_first([chrsz])
        File restriction_sites_ = if defined(restriction_sites) then select_first([restriction_sites]) else select_first([make_restriction_site_locations.restriction_site_locations])
        Int num_bioreps = if defined(bams) then length(select_first([bams])) else length(fastq)

        scatter(i in range(num_bioreps)) {
            # Align fastqs if input bams were not provided
            if (!defined(bams) && defined(reference_index)) {
                scatter(j in range(length(fastq[i]))) {
                    call align { input:
                        fastqs = fastq[i][j],
                        chrsz = chrsz_,
                        idx_tar = select_first([reference_index]),
                        ligation_site = ligation_site,
                        cpu = cpu,
                    }
                    if (length(read_groups) > 0) {
                        call add_read_group_to_bam { input:
                            bam = align.result,
                            read_group = read_groups[i][j],
                        }
                    }
                }
                Array[File] aligned_bams = select_all(if length(read_groups) > 0 then add_read_group_to_bam.bam_with_read_group else align.result)
            }

            Array[File] rep_bam_files = if defined(bams) then select_first([bams])[i] else select_first([aligned_bams])
            Array[File] rep_ligation_counts = if defined(ligation_counts) then select_first([ligation_counts])[i] else select_first([align.norm_res])

            # Scatter across all the bams in the biorep (one per tech rep)
            scatter(j in range(length(rep_bam_files))) {
                call fragment { input:
                    bam_file = rep_bam_files[j],
                    norm_res_input = rep_ligation_counts[j],
                    restriction = restriction_sites_,
                    include_mapq0_reads = include_mapq0_reads,
                }
            }

            Array[File] unmapped = fragment.unmapped
            Array[File] mapq0 = fragment.mapq0
            Array[File] alignable = fragment.alignable

            Array[Array[File]] bams_to_merge = [unmapped, mapq0, alignable]

            scatter(j in range(length(bams_to_merge))){
                call merge { input:
                    bam_files = bams_to_merge[j]
                }
            }

            call merge_sort { input:
                sort_files_ = fragment.sort_file
            }

            # we can collect the alignable.bam using the array merge.out_file
            call dedup { input:
                merged_sort = merge_sort.out_file,
                ligation_site = ligation_site,
                restriction_sites = restriction_sites_,
                alignable_bam = merge.merged_output[2]
            }

            # convert alignable bam to pairs to be consistent with 4DN
            if ( !no_bam2pairs ) {
                call bam2pairs { input:
                    bam_file = dedup.deduped_bam,
                    chrsz_  = chrsz_
                }
            }
        }
    }

    if (defined(input_dedup_pairs) || defined(dedup.out_file)) {
        call merge_pairs_file { input:
            not_merged_pe = select_first([input_dedup_pairs, dedup.out_file])
        }
    }

    if (defined(fragment.alignment_stats) && defined(dedup.library_complexity) && !defined(input_dedup_pairs)) {
        call merge_stats { input:
            alignment_stats = flatten(select_first([fragment.alignment_stats])),
            library_stats = select_first([dedup.library_complexity])
        }
    }

    if (defined(input_dedup_pairs) && defined(alignment_stats) && defined(library_stats)) {
        call merge_stats as merge_stats_from_library_entrypoint { input:
            alignment_stats = select_first([alignment_stats]),
            library_stats = select_first([library_stats]),
        }
    }

    Array[String] qualities = if !defined(input_hic) then DEFAULT_HIC_QUALITIES else []
    if (defined(chrsz) && has_restriction_sites) {
        scatter(i in range(length(qualities))) {
            call create_hic { input:
                pairs_file = select_first([input_pairs, merge_pairs_file.out_file]),
                restriction_sites = if defined(restriction_sites) then select_first([restriction_sites]) else select_first([make_restriction_site_locations.restriction_site_locations]),
                ligation_site = ligation_site,
                chrsz_ = select_first([chrsz]),
                quality = qualities[i],
                assembly_name = assembly_name
            }
        }
    }

    if ( (defined(input_hic) || defined(create_hic.inter)) && !no_call_tads ) {
        call arrowhead { input:
            hic_file = if defined(input_hic) then select_first([input_hic]) else select_first([create_hic.inter])[1]
        }
    }

    if ( (defined(input_hic) || defined(create_hic.inter)) && !no_call_loops ) {
        call hiccups { input:
            hic_file = if defined(input_hic) then select_first([input_hic]) else select_first([create_hic.inter])[1]
        }
    }

    output {
        File? restriction_site_locations = make_restriction_site_locations.restriction_site_locations
        Array[File]? alignable_bam = dedup.deduped_bam
        Array[File?]? out_pairs = bam2pairs.out_file
        Array[File]? out_dedup =  dedup.out_file
        Array[File]? library_complexity_stats_json = dedup.library_complexity_json
        Array[File]? stats = dedup.stats_json
        Array[Array[File]]? alignment_stats_ = fragment.alignment_stats
        Array[Array[File?]?]? bams_with_read_group = add_read_group_to_bam.bam_with_read_group
        File? merged_stats = merge_stats_from_library_entrypoint.merged_stats
        File? merged_stats_json = merge_stats_from_library_entrypoint.merged_stats_json
        File? out_hic_1 = if defined(create_hic.inter) then select_first([create_hic.inter])[0] else input_hic
        File? out_hic_30 = if defined(create_hic.inter) then select_first([create_hic.inter])[1] else input_hic
        File? out_tads = arrowhead.out_file
    }
}

task make_restriction_site_locations {
    input {
        File reference_fasta
        String assembly_name
        String restriction_enzyme
    }

    command <<<
        python3 "$(which generate_site_positions.py)" ~{restriction_enzyme} ~{assembly_name} ~{reference_fasta}
        gzip -n "~{assembly_name}_~{restriction_enzyme}.txt"
    >>>

    output {
        File restriction_site_locations = "~{assembly_name}_~{restriction_enzyme}.txt.gz"
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
        File idx_tar        # reference bwa index tar
        Array[File] fastqs  # [read_end_id]
        File chrsz          # chromosome sizes file
        String ligation_site
        Int? cpu = 32
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
        curr_ostem="result"
        ligation=${"'" + ligation_site + "''"}
        file1=${fastqs[0]}
        file2=${fastqs[1]}
        #count ligations
        source /opt/scripts/common/countligations.sh
        # Align reads
        echo "Running bwa command"
        bwa mem -SP5M -t ${cpu} $reference_index_path ${fastqs[0]} ${fastqs[1]} | samtools view -hbS - > result.bam
    }

    output {
        File result = glob("result.bam")[0]
        File norm_res = glob("result_norm.txt.res.txt")[0]
     }

    runtime {
        cpu : "32"
        memory: "64 GB"
        disks: "local-disk 1000 HDD"
    }
}

task add_read_group_to_bam {
    input {
        File bam
        String read_group
        Int? num_cpus = 8
    }

    command {
        samtools addreplacerg -r "~{read_group}" -o with_read_group.bam --threads "~{num_cpus - 1}" "~{bam}"
    }

    output {
        File bam_with_read_group = "with_read_group.bam"
    }

    runtime {
        cpu : "~{num_cpus}"
        memory: "32 GB"
        disks: "local-disk 1000 HDD"
    }
}

task fragment {
    input {
        File bam_file
        File norm_res_input
        File restriction    # restriction enzyme sites in the reference genome
        Boolean include_mapq0_reads = false
        Int ram_gb = 16
        Float ram_pct = 0.9
    }

    command {
        set -euo pipefail
        samtools view -h ${bam_file} | awk -v "fname"=result -v "mapq0_reads_included"=${if include_mapq0_reads then 1 else 0} -f $(which chimeric_blacklist.awk)

        # if any normal reads were written, find what fragment they correspond
        # to and store that

        echo "Running fragment"
        fragment result_norm.txt result_frag.txt ${restriction}
        echo $(ls)

        # no restriction site !!!!
        # need to add support

        # qc for alignment portion
        cat ${norm_res_input} | awk -f  $(which stats_sub.awk) >> alignment_stats.txt
        paste -d "" ${norm_res_input} > result.res.txt
        python3 $(which jsonify_stats.py) --alignment-stats alignment_stats.txt

        # convert sams to bams and delete the sams
        echo "Converting sam to bam"
        samtools view -hb result_unmapped.sam > unmapped.bam
        rm result_unmapped.sam
        samtools view -hb result_mapq0.sam > mapq0.bam
        rm result_mapq0.sam
        samtools view -hb result_alignable.sam > alignable.bam
        rm result_alignable.sam
        #removed all sam files
        ##restriction used to be site_file

        # sort by chromosome, fragment, strand, and position
        sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S "~{round(ram_pct * ram_gb)}G" result_frag.txt > sort.txt
        gzip -n sort.txt
    }

    output {
        File unmapped = "unmapped.bam"
        File mapq0 = "mapq0.bam"
        File alignable = "alignable.bam"
        File sort_file = "sort.txt.gz"
        File norm_res = "result.res.txt"
        File alignment_stats = "alignment_stats.txt"
        File alignment_stats_json = "alignment_stats.json"
    }

    runtime {
        cpu : "1"
        disks: "local-disk 1000 HDD"
        memory: "~{ram_gb} GB"
    }
}

task merge {
    input {
        Array[File] bam_files
    }

    command {
        set -euo pipefail
        samtools merge merged_bam_files.bam ~{sep=' ' bam_files}
    }

    output {
        File merged_output= "merged_bam_files.bam"
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"

    }
}

task merge_sort {
    input {
        Array[File] sort_files_
        Int ram_gb = 16
        Float ram_pct = 0.9
    }

    command <<<
        set -euo pipefail
        SORT_FILES=sort_files
        mkdir "${SORT_FILES}"
        # Add random number to make filenames unique
        for i in ~{sep=' ' sort_files_}; do mv $i "${SORT_FILES}/sort_${RANDOM}.txt.gz"; done
        # Task test doesn't pass consistently without force option -f due to symlinking
        gzip -df "${SORT_FILES}"/*
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S "~{round(ram_pct * ram_gb)}G" "${SORT_FILES}"/* > merged_sort.txt
        gzip -n merged_sort.txt
    >>>

    output {
        File out_file = "merged_sort.txt.gz"
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
        memory: "~{ram_gb} GB"
    }
}

task dedup {
    input {
        File merged_sort
        String ligation_site
        File restriction_sites
        File alignable_bam
    }

    # Can't use regular {} for command block, parser complains once hits awk command
    command <<<
        set -euo pipefail
        MERGED_SORT_FILE=merged_sort.txt
        touch dups.txt
        touch optdups.txt
        touch merged_nodups.txt
        gzip -dc ~{merged_sort} > $MERGED_SORT_FILE
        awk -f $(which dups.awk) $MERGED_SORT_FILE
        pcr=$(wc -l dups.txt | awk '{print $1}')
        unique=$(wc -l merged_nodups.txt | awk '{print $1}')
        opt=$(wc -l optdups.txt | awk '{print $1}')
        java -jar -Ddevelopment=false /opt/scripts/common/juicer_tools.jar LibraryComplexity $unique $pcr $opt > library_complexity.txt
        /opt/scripts/common/statistics.pl -s ~{restriction_sites} -l ~{"'" + ligation_site + "''"} merged_nodups.txt
        python3 $(which jsonify_stats.py) --library-complexity library_complexity.txt
        python3 $(which jsonify_stats.py) --library-stats stats.txt
        awk '{split($(NF-1), a, "$"); split($NF, b, "$"); print a[3],b[3] > a[2]"_dedup"}' merged_nodups.txt
        samtools view -h ~{alignable_bam} | awk 'BEGIN{OFS="\t"}FNR==NR{for (i=$1; i<=$2; i++){a[i];} next}(!(FNR in a) && $1 !~ /^@/){$2=or($2,1024)}{print}' result_dedup - > result_alignable_dedup.sam
        samtools view -hb result_alignable_dedup.sam > result_alignable_dedup.bam
        rm result_alignable_dedup.sam
        rm ~{alignable_bam}
        gzip -n merged_nodups.txt
    >>>

    output {
        File deduped_bam = "result_alignable_dedup.bam"
        File out_file = "merged_nodups.txt.gz"
        File library_complexity = "library_complexity.txt"
        File library_complexity_json = "library_complexity.json"
        File stats = "stats.txt"
        File stats_json = "stats.json"
        File stats_hists = "stats_hists.m"
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
        memory: "16 GB"
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
        disks: "local-disk 1000 HDD"
    }
}


task merge_pairs_file {
    input {
        Array[File] not_merged_pe
        Int ram_gb = 16
        Float ram_pct = 0.9
    }

    command <<<
        set -euo pipefail
        NOT_MERGED_PE_FILES=not_merged_pes
        mkdir "${NOT_MERGED_PE_FILES}"
        for i in ~{sep=' ' not_merged_pe}; do mv $i "${NOT_MERGED_PE_FILES}/merged_nodups_${RANDOM}.txt.gz"; done
        gzip -df "${NOT_MERGED_PE_FILES}"/*
        sort -m -k2,2d -k6,6d --parallel=8 -S "~{round(ram_pct * ram_gb)}G" "${NOT_MERGED_PE_FILES}"/* > merged_pairs.txt
        gzip -n merged_pairs.txt
    >>>

    output {
        File out_file = "merged_pairs.txt.gz"
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
        memory: "~{ram_gb} GB"
    }
}


task merge_stats {
    # Merge QC statistics from multiple libraries
    input {
        Array[File] alignment_stats
        Array[File] library_stats
    }

    command {
        set -euo pipefail
        awk -f $(which makemega_addstats.awk) ${sep=' ' alignment_stats} ${sep=' ' library_stats} > merged_stats.txt
        python3 $(which jsonify_stats.py) --alignment-stats merged_stats.txt
    }

    output {
        File merged_stats = "merged_stats.txt"
        File merged_stats_json = "merged_stats.json"
    }
}

task create_hic {
    input {
        String ligation_site
        File pairs_file
        File chrsz_
        File restriction_sites
        String quality
        String? assembly_name
    }

    command {
        set -euo pipefail
        MERGED_PAIRS_FILE=merged_pairs.txt
        gzip -dc ~{pairs_file} > $MERGED_PAIRS_FILE
        statistics.pl -q ${quality} -o stats_${quality}.txt -s ${restriction_sites} -l ${"'" + ligation_site + "''"} $MERGED_PAIRS_FILE
        # If the assembly name is empty, then we write chrsz path into file as usual, otherwise, use the assembly name instead of the path
        juicer_tools pre \
            -s stats_${quality}.txt \
            -g stats_${quality}_hists.m \
            -q ${quality} \
            ~{if defined(assembly_name) then "-y " + assembly_name else ""} \
            $MERGED_PAIRS_FILE \
            inter_${quality}.hic \
            ${chrsz_}
        python3 $(which jsonify_stats.py) --alignment-stats stats_${quality}.txt
    }

    output {
        File inter = glob('inter*.hic')[0]
        File stats = glob('stats*.txt')[0]
        File stats_json = glob('stats*.json')[0]
        File stats_hists = glob('stats*hists.m')[0]
    }

    runtime {
        cpu : "1"
        disks: "local-disk 1000 HDD"
        memory : "64 GB"
    }
}


task arrowhead {
    input {
        File hic_file
    }

    command {
        set -euo pipefail
        juicer_tools arrowhead ${hic_file} contact_domains
        gzip -n contact_domains/*
    }

    output {
        File out_file = glob('contact_domains/*.bedpe.gz')[0]
    }

    runtime {
    }
}

task hiccups{
    input {
        File hic_file
    }

    command {
        set -euo pipefail
        java -jar -Ddevelopment=false /opt/scripts/common/juicer_tools.jar hiccups --ignore_sparsity ${hic_file} loops
        gzip -n loops/*.bedpe
    }

    output {
        File out_file = glob("loops/*.bedpe.gz")[0]
    }

    runtime {
        bootDiskSizeGb: "20"
        docker: "encodedcc/hic-pipeline:0.1.0_hiccups"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
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
}
