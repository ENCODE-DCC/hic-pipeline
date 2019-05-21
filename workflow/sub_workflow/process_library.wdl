workflow process_library {
    Array[Array[File]] sub_fastq
    String restriction_enzyme
    File restriction_sites
    File chrsz
    File reference_index
    Int? cpu

    Int fastqs_len = length(sub_fastq)

    # The ligation junctions are consistent with mega.sh
    Map[String, String] restriction_enzyme_to_site = read_map("workflow/restriction_enzyme_to_site.tsv")
    String ligation_site = restriction_enzyme_to_site[restriction_enzyme]

    scatter(j in range(fastqs_len)){
        call align { input:
            fastqs = sub_fastq[j],
            chrsz = chrsz,
            idx_tar = reference_index,
            ligation_site = ligation_site,
            cpu = cpu
        }

        call fragment { input:
            bam_file = align.result,
            norm_res_input = align.norm_res,
            restriction = restriction_sites
        }
    }

    Array[File] collisions = fragment.collisions  #for separate user entry point
    Array[File] collisions_low = fragment.collisions_low_mapq
    Array[File] unmapped = fragment.unmapped
    Array[File] mapq0 = fragment.mapq0
    Array[File] alignable = fragment.alignable

    Array[Array[File]] bams_to_merge = [collisions, collisions_low, unmapped, mapq0, alignable]

    scatter(i in range(5)){
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
        restriction_sites = restriction_sites,
        alignable_bam = merge.merged_output[4]
    }
    
    # convert alignable bam to pairs to be consistent with 4DN
    call bam2pairs { input:
        bam_file = dedup.deduped_bam,
        chrsz_  = chrsz
    }

    output {
        File alignable_bam = dedup.deduped_bam
        File pairs_file = bam2pairs.out_file
        File library_dedup = dedup.out_file
        File stats_json = dedup.stats_json
        File library_stats = dedup.library_complexity
        File library_stats_json = dedup.library_complexity_json
        File stats = dedup.stats
        File stats_hists = dedup.stats_hists
        Array[File] alignments_stats = fragment.stats_sub_result
        Array[File] alignment_stats_json = fragment.stats_sub_result_json
    }
}

task align {
    File idx_tar        # reference bwa index tar
    Array[File] fastqs  # [read_end_id]
    File chrsz          # chromosome sizes file
    String ligation_site
    Int? cpu

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
        bwa mem -SP5M -t ${select_first([cpu,32])} $reference_index_path ${fastqs[0]} ${fastqs[1]} | samtools view -hbS - > result.bam
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
    File bam_file
    File norm_res_input
    File restriction    # restriction enzyme sites in the reference genome

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
            exit 1x
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
        File stats_sub_result = glob("alignment_stats.txt")[0]
        File stats_sub_result_json = glob("alignment_stats.json")[0]
    }

    runtime {
        cpu : "1"
        disks: "local-disk 1000 HDD"
    }
}

task merge {
    Array[File] bam_files
    
    command <<<
        samtools merge merged_bam_files.bam ${sep=' ' bam_files} 
    >>>

    output {
        File merged_output= glob('merged_bam_files.bam')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
        
    }
}

task merge_sort {
   Array[File] sort_files_

    command {
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 90% ${sep=' ' sort_files_}  > merged_sort.txt
    }

    output {
        File out_file = glob('merged_sort.txt')[0]
    }

    runtime {
        cpu : "8"
        disks: "local-disk 1000 HDD"
    }
}

task dedup {
    File merged_sort
    String ligation_site
    File restriction_sites
    File alignable_bam

    command <<<
        touch dups.txt
        touch optdups.txt
        touch merged_nodups.txt
        awk -f /opt/scripts/common/dups.awk ${merged_sort}
        pcr=$(wc -l dups.txt | awk '{print $1}')
        unique=$(wc -l merged_nodups.txt | awk '{print $1}')
        opt=$(wc -l optdups.txt | awk '{print $1}')
        java -jar -Ddevelopment=false /opt/scripts/common/juicer_tools.jar LibraryComplexity $unique $pcr $opt > library_complexity.txt
        /opt/scripts/common/statistics.pl -s ${restriction_sites} -l ${ligation_site} merged_nodups.txt
        python3 /opt/hic-pipeline/src/jsonify_stats.py --library-complexity library_complexity.txt
        python3 /opt/hic-pipeline/src/jsonify_stats.py --library-stats stats.txt
        awk '{split($(NF-1), a, "$"); split($NF, b, "$"); print a[3],b[3] > a[2]"_dedup"}' merged_nodups.txt
        samtools view -h ${alignable_bam} | awk 'BEGIN{OFS="\t"}FNR==NR{for (i=$1; i<=$2; i++){a[i];} next}(!(FNR in a) && $1 !~ /^@/){$2=or($2,1024)}{print}' result_dedup - > result_alignable_dedup.sam
        samtools view -hb result_alignable_dedup.sam > result_alignable_dedup.bam
        rm result_alignable_dedup.sam
        rm ${alignable_bam}
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
    File bam_file
    File chrsz_

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
