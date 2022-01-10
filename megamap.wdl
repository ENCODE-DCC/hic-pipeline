version 1.0

import "./hic.wdl"

workflow megamap {
    meta {
        version: "1.8.0"
        caper_docker: "encodedcc/hic-pipeline:1.8.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.8.0"
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

    output {
        File fasta_index = create_reference_fasta_index.fasta_index
        File sequence_dictionary = create_gatk_references.sequence_dictionary
        File interval_list = create_gatk_references.interval_list
        File snp_vcf = gatk.snp_vcf
        File indel_vcf = gatk.indel_vcf
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
