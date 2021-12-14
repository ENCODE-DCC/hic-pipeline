version 1.0

import "./hic.wdl"

workflow megamap {
    meta {
        version: "1.10.0"
        caper_docker: "encodedcc/hic-pipeline:1.10.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.10.0"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json"
    }

    input {
        File reference_fasta
        Array[File] bams
        # From GATK bundle
        # https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
        # https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
        File dbsnp_vcf
        File dbsnp_vcf_index
        File hapmap_vcf_index
        File hapmap_vcf
        File mills_vcf
        File mills_vcf_index
        File omni_vcf
        File omni_vcf_index
        Int? gatk_num_cpus
        # Only for testing purposes
        Boolean no_bundle = false
    }

    call hic.merge as merged { input:
        bams = bams,
    }

    call hic.bam_to_pre as converted { input:
        bam = merged.bam,
        quality = 1,
        # merged_nodups = True?,
    }

    call create_fasta_index as create_reference_fasta_index { input:
        fasta = reference_fasta,
    }

    call create_gatk_references { input:
        reference_fasta = reference_fasta,
        reference_fasta_index = create_reference_fasta_index.fasta_index,
        output_stem = basename(reference_fasta, ".fasta.gz")
    }

    call gatk { input:
        bam = merged.bam,
        reference_fasta = reference_fasta,
        reference_fasta_index = create_reference_fasta_index.fasta_index,
        sequence_dictionary = create_gatk_references.sequence_dictionary,
        interval_list = create_gatk_references.interval_list,
        mills_vcf = mills_vcf,
        omni_vcf = omni_vcf,
        hapmap_vcf = hapmap_vcf,
        dbsnp_vcf = dbsnp_vcf,
        mills_vcf_index = mills_vcf_index,
        omni_vcf_index = omni_vcf_index,
        hapmap_vcf_index = hapmap_vcf_index,
        dbsnp_vcf_index = dbsnp_vcf_index,
        no_bundle = no_bundle,
        num_cpus = gatk_num_cpus,
    }

    call run_3d_dna { input:
        reference_fasta = reference_fasta,
        merged_nodups = converted.pre,
    }
}

task create_fasta_index {
    input {
        File fasta
    }

    command <<<
        set -euo pipefail
        gzip -dc ~{fasta} > ~{basename(fasta, ".gz")}
        samtools faidx ~{basename(fasta, ".gz")}
    >>>

    output {
        File fasta_index = "~{basename(fasta, '.gz')}.fai"
    }
}

task create_gatk_references {
    input {
        File reference_fasta
        File reference_fasta_index
        String output_stem
    }

    command <<<
        set -euo pipefail
        gzip -dc ~{reference_fasta} > ~{basename(reference_fasta, ".gz")}
        mv ~{reference_fasta_index} .
        gatk \
            CreateSequenceDictionary \
            --REFERENCE ~{basename(reference_fasta, ".gz")} \
            --OUTPUT ~{output_stem}.dict \
            --URI ~{basename(reference_fasta, ".gz")}
        gatk \
            ScatterIntervalsByNs \
            --REFERENCE ~{basename(reference_fasta, ".gz")} \
            --OUTPUT_TYPE ACGT \
            --MAX_TO_MERGE 500 \
            --OUTPUT ~{output_stem}.interval_list
    >>>

    output {
        File sequence_dictionary = "~{output_stem}.dict"
        File interval_list = "~{output_stem}.interval_list"
    }

    runtime {
        cpu : "1"
        memory: "16 GB"
        disks: "local-disk 100 HDD"
    }
}

task gatk {
    input {
        File bam
        File reference_fasta
        File reference_fasta_index
        File sequence_dictionary
        File interval_list
        File mills_vcf
        File omni_vcf
        File hapmap_vcf
        File dbsnp_vcf
        File mills_vcf_index
        File omni_vcf_index
        File hapmap_vcf_index
        File dbsnp_vcf_index
        # Only for testing purposes
        Boolean no_bundle = false
        Int num_cpus = 16
    }

    String final_snp_vcf_name = "snp.out.vcf"
    String final_indel_vcf_name = "indel.out.vcf"

    command <<<
        mkdir bundle
        if [[ ~{if(no_bundle) then "0" else "1"} -eq 1 ]]
        then
            mv \
                ~{mills_vcf} \
                ~{omni_vcf} \
                ~{hapmap_vcf} \
                ~{dbsnp_vcf} \
                ~{mills_vcf_index} \
                ~{omni_vcf_index} \
                ~{hapmap_vcf_index} \
                ~{dbsnp_vcf_index} \
                bundle
        fi
        mkdir reference
        mv ~{reference_fasta_index} ~{sequence_dictionary} ~{interval_list} reference
        gzip -dc ~{reference_fasta} > reference/~{basename(reference_fasta, ".gz")}
        run-gatk-after-juicer2.sh \
            -r reference/~{basename(reference_fasta, ".gz")} \
            ~{if !no_bundle then "--gatk-bundle bundle" else ""} \
            --threads ~{num_cpus} \
            ~{bam}
        gzip -n ~{final_snp_vcf_name}
        gzip -n ~{final_indel_vcf_name}
    >>>

    output {
        File snp_vcf = "~{final_snp_vcf_name}.gz"
        File indel_vcf = "~{final_indel_vcf_name}.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        memory: "128 GB"
        disks: "local-disk 1000 HDD"
    }
}


task run_3d_dna {
    input {
        File reference_fasta
        File merged_nodups
    }

    command <<<
        set -euo pipefail
        export MERGED_NODUPS_FILENAME=merged_nodups.txt
        gzip -dc ~{merged_nodups} > $MERGED_NODUPS_FILENAME
        bash /opt/3d-dna/run-3ddna-pipeline.sh ~{reference_fasta} $MERGED_NODUPS_FILENAME
        ls
    >>>

    output {
        # fasta files
        # chromosome-length scaffolds, small and tiny scaffolds
        File scaffolds_fasta = "HiC.fasta"

        # .hic files
        # sandboxed contact map corresponding to the HiC.fasta reference
        File hic = "HiC.hic"
        # after sealing stage
        File after_sealing_rawchrom_hic = "rawchrom.hic"
        File after_sealing_final_hic = "final.hic"
        # after polishing stage
        File polished_hic = "polished.hic"
        # after editing and scaffolding
        File resolved_hic = "resolved.hic"

        # Scaffold boundary files (Juicebox 2D annotation format)
        File scaffold_track = "scaffold_track.txt"
        File superscaffold_track = "superscaf_track.txt"

        # Tracks illustrating putative misjoins;
        File bed = "FINAL.bed"
        File wig = "FINAL.wig"

        # .assembly (supersedes .cprops and .asm files)
        # Custom file format that tracks modifications to the input contigs at various stages in the assembly
        # after the addition of gaps to the chromosome-length assembly;
        File assembly = "HiC.assembly"
        # after sealing stage
        File after_sealing_final_assembly = "final.assembly"
        # after polishing stage
        File polished_assembly = "polished.assembly"
        # after editing and scaffolding
        File resolved_assembly = "resolved.assembly"

        # supplementary files
        # list of problematic regions (Juicebox 2D annotation format)
        Array[File] edits_for_step = glob("edits.for.step.*.txt")
        Array[File] mismatches_at_step = glob("mismatches.at.step.*.txt")
        Array[File] suspect_2D_at_step = glob("suspect_2D.at.step.*.txt")
        # pairwise alignment data for alternative haplotype candidates, parsed from LASTZ output.
        File pairwise_alignments = "alignments.txt"
    }

    runtime {
        cpu : "16"
        disks: "local-disk 200 HDD"
        memory: "32 GB"
    }
}
