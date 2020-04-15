version 1.0

#CAPER docker quay.io/encode-dcc/hic-pipeline:0.1.0
#CAPER singularity docker://quay.io/encode-dcc/hic-pipeline:0.1.0

workflow hic {
    input {
        # Main entrypoint, need to specify all five of these values when running from fastqs
        Array[Array[Array[File]]] fastq = []
        String restriction_enzyme
        File? restriction_sites
        File? chrsz
        File? reference_index

        # Entrypoint right before hic generation
        File? input_pairs

        # Entrypoint for loop and TAD calls
        File? input_hic

        # Inputs and logic for entrypoint after library processing
        Array[File]? input_dedup_pairs
        Array[File?] library_stats = []
        Array[File?] library_stats_hists = []

        Boolean? no_bam2pairs = false
        Boolean? no_call_loops = false
        Boolean? no_call_tads = false
        Boolean? include_mapq0_reads = false
        Array[String]? input_ligation_junctions
        Int cpu = 32
        String? assembly_name
    }

    parameter_meta {
        fastq: "Twice nested array of input fastqs, takes form of [lib_id][fastq_id][read_end_id]"
        restriction_enzyme: "The name of the restriction enzyme used to generate the Hi-C libraries"
        restriction_sites: "A text file containing cut sites for the given restriction enzyme. You should generate this file using this script: https://github.com/aidenlab/juicer/blob/encode/misc/generate_site_positions.py"
        chrsz: "A chromosome sizes file for the desired assembly, this is a tab-separated text file whose rows take the form [chromosome] [size]"
        reference_index: "A pregenerated BWA index for the desired assembly"
        input_pairs: "A text file containing the paired fragments to use to generate the .hic contact maps, a detailed format description can be found here: https://github.com/aidenlab/juicer/wiki/Pre#long-format"
        input_hic: "An input .hic file for which to call loops and domains"
        input_dedup_pairs: "An array consisting of text files of paired fragments, one per library, same format as input_pairs"
        library_stats: "An array of library statistics text files, one per library"
        library_stats_hists: "An array of library statistics .m files, one per library"
        input_ligation_junctions: "An array of ligation sites, useful for megamaps"
        no_bam2pairs: "If set to `true`, avoid generating .pairs files, defaults to false"
        no_call_loops: "If set to `true`, avoid calling loops with hiccups, defaults to false"
        no_call_tads: "If set to `true`, avoid calling domains with arrowhead, defaults to false"
        include_mapq0_reads: "If set to `true`, chimeric reads (3 alignments) with one MAPQ 0 read will be classified as normal paired reads with the MAPQ 0 read discarded. If `false`, such reads will be classified as low mapq collisions, defaults to false"
        cpu: "Number of threads to use for bwa alignment"
        assembly_name: "Name of assembly to insert into hic file header, recommended to specify for reproducbility otherwise hic file will be nondeterministic"
    }

    # Pipeline internal "global" variables: do not specify as input
    # These ligation junctions are consistent with mega.sh
    Map[String, String] RESTRICTION_ENZYME_TO_SITE = {
        "HindIII": "AAGCTAGCTT",
        "DpnII": "GATCGATC",
        "MboI": "GATCGATC",
    }

    # Default MAPQ thresholds for generating .hic contact maps
    Array[String] DEFAULT_HIC_QUALITIES = ["1", "30"]

    # Determine ligation site from enzyme name
    String ligation_site = RESTRICTION_ENZYME_TO_SITE[restriction_enzyme]
    # Prepare array of restriction sites for megamap
    Array[String] ligation_junctions = select_first([input_ligation_junctions, [ligation_site]])

    # scatter over libraries
    if (defined(chrsz) && defined(reference_index) && defined(restriction_sites)) {
        # A WDL technicality here.
        File chrsz_ = select_first([chrsz])
        File reference_index_ = select_first([reference_index])

        scatter(i in range(length(fastq))) {
            scatter(j in range(length(fastq[i]))) {
                call align { input:
                    fastqs = fastq[i][j],
                    chrsz = chrsz_,
                    idx_tar = reference_index_,
                    ligation_site = ligation_site,
                    cpu = cpu,
                }

                call fragment { input:
                    bam_file = align.result,
                    norm_res_input = align.norm_res,
                    restriction = select_first([restriction_sites])
                }
            }

            Array[File] collisions = fragment.collisions  #for separate user entry point
            Array[File] collisions_low = fragment.collisions_low_mapq
            Array[File] unmapped = fragment.unmapped
            Array[File] mapq0 = fragment.mapq0
            Array[File] alignable = fragment.alignable

            Array[Array[File]] bams_to_merge = [collisions, collisions_low, unmapped, mapq0, alignable]

            scatter(i in range(length(bams_to_merge))){
                call merge { input:
                    bam_files = bams_to_merge[i]
                }
            }

            call merge_sort { input:
                sort_files_ = fragment.sort_file
            }

            # we can collect the alignable.bam using the array merge.out_file
            call dedup { input:
                merged_sort = merge_sort.out_file,
                ligation_site = ligation_site,
                restriction_sites = select_first([restriction_sites]),
                alignable_bam = merge.merged_output[4]
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

    if (defined(fragment.alignment_stats) && defined(dedup.library_complexity)) {
        call merge_stats { input:
            alignment_stats = flatten(select_first([fragment.alignment_stats])),
            library_stats = select_first([dedup.library_complexity])
        }
    }

    Array[String] qualities = if !defined(input_hic) then DEFAULT_HIC_QUALITIES else []
    if (defined(chrsz) && defined(restriction_sites)) {
        scatter(i in range(length(qualities))) {
            call create_hic { input:
                pairs_file = select_first([input_pairs, merge_pairs_file.out_file]),
                restriction_sites = select_first([restriction_sites]),
                ligation_junctions = ligation_junctions,
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
        Array[File]? alignable_bam = dedup.deduped_bam
        Array[File?]? out_pairs = bam2pairs.out_file
        Array[File]? out_dedup =  dedup.out_file
        Array[File]? library_complexity_stats_json = dedup.library_complexity_json
        Array[File]? stats = dedup.stats_json
        Array[Array[File]]? alignment_stats = fragment.alignment_stats
        Array[Array[File]]? alignment_stats_json = fragment.alignment_stats_json
        File? out_hic_1 = if defined(create_hic.inter) then select_first([create_hic.inter])[0] else input_hic
        File? out_hic_30 = if defined(create_hic.inter) then select_first([create_hic.inter])[1] else input_hic
        File? out_tads = arrowhead.out_file
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
        ligation=${ligation_site}
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

task fragment {
    input {
        File bam_file
        File norm_res_input
        File restriction    # restriction enzyme sites in the reference genome
        Boolean include_mapq0_reads = false
    }

    command {
        samtools view -h ${bam_file} | awk -v "fname"=result -v "mapq0_reads_included"=${if include_mapq0_reads then 1 else 0} -f $(which chimeric_blacklist.awk)

        # if any normal reads were written, find what fragment they correspond
        # to and store that

        echo "Running fragment"
        fragment result_norm.txt result_frag.txt ${restriction}
        echo $(ls)

        # no restriction site !!!!
        # need to add support

        # qc for alignment portion
        cat ${norm_res_input} *.res.txt | awk -f  $(which stats_sub.awk) >> alignment_stats.txt
        paste -d "" ${norm_res_input} *.res.txt > result.res.txt
        python3 $(which jsonify_stats.py) --alignment-stats alignment_stats.txt

        # convert sams to bams and delete the sams
        echo "Converting sam to bam"
        samtools view -hb result_collisions.sam > collisions.bam
        rm result_collisions.sam
        samtools view -hb result_collisions_low_mapq.sam > collisions_low_mapq.bam
        rm result_collisions_low_mapq.sam
        samtools view -hb result_unmapped.sam > unmapped.bam
        rm result_unmapped.sam
        samtools view -hb result_mapq0.sam > mapq0.bam
        rm result_mapq0.sam
        samtools view -hb result_alignable.sam > alignable.bam
        rm result_alignable.sam
        #removed all sam files
        ##restriction used to be site_file

        # sort by chromosome, fragment, strand, and position
        sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 90% result_frag.txt > sort.txt
        if [ $? -ne 0 ]
        then
            echo "***! Failure during sort"
            exit 1
        fi
        gzip -n sort.txt
    }

    output {
        File collisions = glob("collisions.bam")[0]
        File collisions_low_mapq = glob("collisions_low_mapq.bam")[0]
        File unmapped = glob("unmapped.bam")[0]
        File mapq0 = glob("mapq0.bam")[0]
        File alignable = glob("alignable.bam")[0]
        File sort_file = glob("sort.txt.gz")[0]
        File norm_res = glob("result.res.txt")[0]
        File alignment_stats = glob("alignment_stats.txt")[0]
        File alignment_stats_json = glob("alignment_stats.json")[0]
    }

    runtime {
        cpu : "1"
        disks: "local-disk 1000 HDD"
        memory: "16 GB"
    }
}

task merge {
    input {
        Array[File] bam_files
    }

    command {
        samtools merge merged_bam_files.bam ~{sep=' ' bam_files}
    }

    output {
        File merged_output= glob('merged_bam_files.bam')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"

    }
}

task merge_sort {
    input {
        Array[File] sort_files_
    }

    command {
        SORT_FILES=sort_files
        mkdir $SORT_FILES
        mv ~{sep=' ' sort_files_} $SORT_FILES
        for i in $SORT_FILES/*; do gzip -d $i; done
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 90% $SORT_FILES/*  > merged_sort.txt
        gzip -n merged_sort.txt
    }

    output {
        File out_file = glob('merged_sort.txt.gz')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
        memory: "16 GB"
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
        /opt/scripts/common/statistics.pl -s ~{restriction_sites} -l ~{ligation_site} merged_nodups.txt
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
        File deduped_bam = glob('result_alignable_dedup.bam')[0]
        File out_file = glob('merged_nodups.txt.gz')[0]
        File library_complexity = glob('library_complexity.txt')[0]
        File library_complexity_json = glob('library_complexity.json')[0]
        File stats = glob('stats.txt')[0]
        File stats_json = glob('stats.json')[0]
        File stats_hists = glob('stats_hists.m')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}

task bam2pairs {
    input {
        File bam_file
        File chrsz_
    }

    command {
        bam2pairs -c ${chrsz_} ${bam_file} pairix
    }

    output {
        File out_file = glob('pairix.*.pairs.gz')[0]
        File pairs_index = glob('pairix.*.pairs.gz.px2')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}


task merge_pairs_file {
    input {
        Array[File] not_merged_pe
    }

    command {
        NOT_MERGED_PE_FILES=not_merged_pes
        mkdir $NOT_MERGED_PE_FILES
        mv ~{sep=' ' not_merged_pe} $NOT_MERGED_PE_FILES
        for i in $NOT_MERGED_PE_FILES/*; do gzip -d $i; done
        sort -m -k2,2d -k6,6d --parallel=8 -S 10% $NOT_MERGED_PE_FILES/* > merged_pairs.txt
        gzip -n merged_pairs.txt
    }

    output {
        File out_file = glob('merged_pairs.txt.gz')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}


task merge_stats {
    # Merge QC statistics from multiple libraries
    input {
        Array[File] alignment_stats
        Array[File] library_stats
    }

    command {
        awk -f $(which makemega_addstats.awk) ${sep=' ' alignment_stats} ${sep=' ' library_stats} > merged_stats.txt
        python3 $(which jsonify_stats.py) --alignment-stats merged_stats.txt
    }

    output {
        File merged_stats = glob('merged_stats.txt')[0]
        File merged_stats_json = glob("merged_stats.json")[0]
    }
}

task create_hic {
    input {
        Array[String] ligation_junctions
        File pairs_file
        File chrsz_
        File restriction_sites
        String quality
        String? assembly_name
    }

    command {
        MERGED_PAIRS_FILE=merged_pairs.txt
        gzip -dc ~{pairs_file} > $MERGED_PAIRS_FILE
        statistics.pl -q ${quality} -o stats_${quality}.txt -s ${restriction_sites} -l ${sep=' ' ligation_junctions} $MERGED_PAIRS_FILE
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
        java -jar -Ddevelopment=false /opt/scripts/common/juicer_tools.jar hiccups --ignore_sparsity ${hic_file} loops
        gzip -n loops/*.bedpe
    }

    output {
        File out_file = glob("loops/*.bedpe.gz")[0]
    }

    runtime {
        bootDiskSizeGb: "20"
        docker: "quay.io/encode-dcc/hiccups:0.1.0"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
    }
}
