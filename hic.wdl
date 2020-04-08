version 1.0

#CAPER docker quay.io/encode-dcc/hic-pipeline:template
#CAPER singularity docker://quay.io/encode-dcc/hic-pipeline:template

workflow hic {
    input {
        Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
        Array[Array[Array[File]]] input_bams = [] #[lib_id[[collisions1,collisions2],[collisions_low],[unmapped],[mapq0],[alignable]], 
        Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
        Array[File] input_merged_sort = []
        File? input_pairs
        File? input_hic
        File? sub_ms
        String? assembly_name

        # Inputs and logic for entrypoint after library processing
        Array[File]? input_dedup_pairs
        Array[File?] library_stats = []
        Array[File?] library_stats_hists = []
        Array[String]? input_ligation_junctions

        # Inputs for library processing
        String restriction_enzyme
        File? restriction_sites
        File? chrsz
        File? reference_index
        Int cpu = 32
        Boolean? no_call_loops = false
        Boolean? no_call_tads = false
        Boolean? no_bam2pairs = false
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

    # Default MAPQ thresholds for .hic contact map generation

    #determine range of scatter
    Int lib_length = if length(fastq) > 0 then length(fastq)
    else if length(input_bams) > 0 then length(input_bams) ##technically the number should be same for bams and sort_files
    else if length(input_sort_files) > 0 then length(input_sort_files)
    else length(input_merged_sort)

    # scatter over libraries
    if (defined(chrsz) && defined(reference_index) && defined(restriction_sites)) {
        # A WDL technicality here.
        File chrsz_ = select_first([chrsz])
        File reference_index_ = select_first([reference_index])

        scatter(i in range(lib_length)) {
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

    File hic_file = select_first([input_hic, select_first([create_hic.inter])[1]])

    if ( (defined(input_hic) || defined(create_hic.inter)) && !no_call_tads ) {
        call arrowhead { input:
            hic_file = hic_file
        }
    }

    if ( (defined(input_hic) || defined(create_hic.inter)) && !no_call_loops ) {
        call hiccups { input:
             hic_file = hic_file
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
    }

    command {
        samtools view -h ${bam_file} | awk -v "fname"=result -f /opt/scripts/common/chimeric_blacklist.awk
 
        # if any normal reads were written, find what fragment they correspond
        # to and store that
        
        echo "Running fragment"
        /opt/scripts/common/fragment.pl result_norm.txt result_frag.txt ${restriction}   
        echo $(ls)

        # no restriction site !!!!
        # need to add support
       
        # qc for alignment portion
        cat ${norm_res_input} *.res.txt | awk -f /opt/scripts/common/stats_sub.awk >> alignment_stats.txt
        paste -d "" ${norm_res_input} *.res.txt > result.res.txt
        python3 /opt/hic-pipeline/src/jsonify_stats.py --alignment-stats alignment_stats.txt

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
    }

    output {
        File collisions = glob("collisions.bam")[0]
        File collisions_low_mapq = glob("collisions_low_mapq.bam")[0]
        File unmapped = glob("unmapped.bam")[0]
        File mapq0 = glob("mapq0.bam")[0]
        File alignable = glob("alignable.bam")[0]
        File sort_file = glob("sort.txt")[0]
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
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 90% ${sep=' ' sort_files_}  > merged_sort.txt
    }

    output {
        File out_file = glob('merged_sort.txt')[0]
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
        touch dups.txt
        touch optdups.txt
        touch merged_nodups.txt
        awk -f /opt/scripts/common/dups.awk ~{merged_sort}
        pcr=$(wc -l dups.txt | awk '{print $1}')
        unique=$(wc -l merged_nodups.txt | awk '{print $1}')
        opt=$(wc -l optdups.txt | awk '{print $1}')
        java -jar -Ddevelopment=false /opt/scripts/common/juicer_tools.jar LibraryComplexity $unique $pcr $opt > library_complexity.txt
        /opt/scripts/common/statistics.pl -s ~{restriction_sites} -l ~{ligation_site} merged_nodups.txt
        python3 /opt/hic-pipeline/src/jsonify_stats.py --library-complexity library_complexity.txt
        python3 /opt/hic-pipeline/src/jsonify_stats.py --library-stats stats.txt
        awk '{split($(NF-1), a, "$"); split($NF, b, "$"); print a[3],b[3] > a[2]"_dedup"}' merged_nodups.txt
        samtools view -h ~{alignable_bam} | awk 'BEGIN{OFS="\t"}FNR==NR{for (i=$1; i<=$2; i++){a[i];} next}(!(FNR in a) && $1 !~ /^@/){$2=or($2,1024)}{print}' result_dedup - > result_alignable_dedup.sam
        samtools view -hb result_alignable_dedup.sam > result_alignable_dedup.bam
        rm result_alignable_dedup.sam
        rm ~{alignable_bam}
    >>>

    output {
        File deduped_bam = glob('result_alignable_dedup.bam')[0]
        File out_file = glob('merged_nodups.txt')[0]
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
        /opt/pairix-0.3.6/util/bam2pairs/bam2pairs -c ${chrsz_} ${bam_file} pairix
    }

    output {
        File out_file = glob('pairix.*.pairs')[0]
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
        sort -m -k2,2d -k6,6d --parallel=8 -S 10% ${sep=' ' not_merged_pe}  > merged_pairs.txt
    }
    
    output {
        File out_file = glob('merged_pairs.txt')[0]
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
        awk -f /opt/scripts/common/makemega_addstats.awk ${sep=' ' alignment_stats} ${sep=' ' library_stats} > merged_stats.txt
        python3 /opt/hic-pipeline/src/jsonify_stats.py --alignment-stats merged_stats.txt
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
        /opt/scripts/common/statistics.pl -q ${quality} -o stats_${quality}.txt -s ${restriction_sites} -l ${sep=' ' ligation_junctions} ${pairs_file}
        ASSEMBLY_NAME=${default='' assembly_name}
        # If the assembly name is empty, then we write chrsz path into file as usual, otherwise, use the assembly name instead of the path
        if [ -z "$ASSEMBLY_NAME" ]; then
            /opt/scripts/common/juicer_tools pre -s stats_${quality}.txt -g stats_${quality}_hists.m -q ${quality} ${pairs_file} inter_${quality}.hic ${chrsz_}
        else
            /opt/scripts/common/juicer_tools pre -s stats_${quality}.txt -g stats_${quality}_hists.m -q ${quality} -y $ASSEMBLY_NAME ${pairs_file} inter_${quality}.hic ${chrsz_}
        fi
        python3 /opt/hic-pipeline/src/jsonify_stats.py --alignment-stats stats_${quality}.txt
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
        /opt/scripts/common/juicer_tools arrowhead ${hic_file} contact_domains
    }

    output {
        File out_file = glob('contact_domains/*.bedpe')[0]
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
    }
    
    output {
        File out_file = glob("loops/*.bedpe")[0]
    }
    
    runtime {
        bootDiskSizeGb: "20"
        docker: "quay.io/encode-dcc/hiccups:template"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
    }
}
