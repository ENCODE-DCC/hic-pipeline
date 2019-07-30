##Encode DCC Hi-C pipeline

#CAPER docker quay.io/encode-dcc/hic-pipeline:template
#CAPER singularity docker://quay.io/encode-dcc/hic-pipeline:template

import "../../workflow/sub_workflow/process_library.wdl" as sub

workflow hic {
    #User inputs 
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[Array[File]]] input_bams = [] #[lib_id[[collisions1,collisions2],[collisions_low],[unmapped],[mapq0],[alignable]], 
    Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
    Array[File] input_merged_sort = []
    File? input_pairs
    File? input_hic
    File? sub_ms
    String? assembly_name

    # Inputs and logic for entrypoint after library processing
    Array[File?] input_dedup_pairs = []
    Array[File?] library_stats = []
    Array[File?] library_stats_hists = []
    Array[String?] input_ligation_junctions = []

    # Inputs for library processing
    String restriction_enzyme
    File? restriction_sites
    File? chrsz
    File? reference_index
    Int? cpu
    Boolean? no_call_loops = false
    Boolean? no_call_tads = false

    #determine range of scatter
    Int lib_length = if length(fastq) > 0 then length(fastq)
    else if length(input_bams) > 0 then length(input_bams) ##technically the number should be same for bams and sort_files
    else if length(input_sort_files) > 0 then length(input_sort_files)
    else length(input_merged_sort)

    # scatter over libraries
    scatter(i in range(lib_length)) {
        call sub.process_library { input:
            sub_fastq = fastq[i],
            chrsz = chrsz,
            reference_index = reference_index,
            restriction_enzyme = restriction_enzyme,
            cpu = cpu,
            restriction_sites = restriction_sites
        }
    }

    call merge_pairs_file { input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else process_library.library_dedup
    }

    call merge_stats { input:
        alignment_stats = flatten(process_library.alignments_stats),
        library_stats = process_library.library_stats
    }

    # Prepare array of restriction sites for megamap
    Map[String, String] restriction_enzyme_to_site = read_map("workflow/restriction_enzyme_to_site.tsv")
    String? ligation_junction = if defined(restriction_enzyme) then restriction_enzyme_to_site[restriction_enzyme] else ""
    Array[String] ligation_junctions = select_all(if defined(restriction_enzyme) then [ligation_junction] else input_ligation_junctions)

    Array[String] qualities = if !defined(input_hic) then ["1", "30"] else []
    scatter(i in range(length(qualities))) {
        call create_hic { input:
            pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
            restriction_sites = restriction_sites,
            ligation_junctions = ligation_junctions,
            chrsz_ = chrsz,
            quality = qualities[i],
            assembly_name = assembly_name
        }
    }

    if ( (defined(input_hic) || defined(create_hic.inter)) && !no_call_tads ) {
        call arrowhead { input:
            hic_file = if defined(input_hic) then input_hic else create_hic.inter[1]
        }
    }

    if ( !no_call_loops ) {
        call hiccups { input:
            hic_file = if defined(input_hic) then input_hic else create_hic.inter[1]
        }
    }

    output {
        # Sub-workflow processing a library outputs
        Array[File] out_merged_align = process_library.alignable_bam
        Array[File] out_pairs = process_library.pairs_file
        Array[File] out_dedup = process_library.library_dedup

        # Create hic outputs
        File? out_hic_1 = if length(create_hic.inter) > 1 then create_hic.inter[0] else input_hic
        File? out_hic_30 = if length(create_hic.inter) > 1 then create_hic.inter[1] else input_hic

        # TADs output
        File? out_tads = arrowhead.out_file
        # HiCCUps output
        File? out_hiccups = hiccups.out_file

        #QC outputs
        Array[File] library_complexity = process_library.library_stats_json
        Array[File] stats = process_library.stats_json
        Array[Array[File]] alignments_stats = process_library.alignments_stats
    }
}


task merge_pairs_file{
    Array[File?] not_merged_pe

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
    Array[File?] alignment_stats
    Array[File?] library_stats

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
    Array[String] ligation_junctions
    File pairs_file
    File chrsz_
    File restriction_sites
    String quality
    String? assembly_name

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
    File hic_file

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
    File hic_file
    
    command {
        java -jar -Ddevelopment=false /opt/scripts/common/juicer_tools.jar hiccups --ignore_sparsity ${hic_file} loops
    }
    
    output {
        File out_file = glob("loops/*.bedpe")[0]
    }
    
    runtime {
        docker: "quay.io/encode-dcc/hiccups:template"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: ["us-east1-b"]
    }     
}
