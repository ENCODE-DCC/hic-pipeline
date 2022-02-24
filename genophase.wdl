version 1.0

import "./hic.wdl"

workflow genophase {
    meta {
        version: "1.11.2"
        caper_docker: "encodedcc/hic-pipeline:1.11.2"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.11.2"
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
        Int? gatk_disk_size_gb
        Int? gatk_ram_gb
        Int? run_3d_dna_num_cpus
        Int? run_3d_dna_disk_size_gb
        Int? run_3d_dna_ram_gb
        Boolean no_phasing = false
        # Only for testing purposes
        Boolean no_bundle = false

        String docker = "encodedcc/hic-pipeline:1.11.2"
        String singularity = "docker://encodedcc/hic-pipeline:1.11.2"
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }

    call hic.merge as merged { input:
        bams = bams,
        runtime_environment = runtime_environment,
    }

    call create_fasta_index as create_reference_fasta_index { input:
        fasta = reference_fasta,
        runtime_environment = runtime_environment,
    }

    call create_gatk_references { input:
        reference_fasta = reference_fasta,
        reference_fasta_index = create_reference_fasta_index.fasta_index,
        output_stem = basename(reference_fasta, ".fasta.gz"),
        runtime_environment = runtime_environment,
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
        ram_gb = gatk_ram_gb,
        disk_size_gb = gatk_disk_size_gb,
        runtime_environment = runtime_environment,
    }

    if (!no_phasing) {
        call run_3d_dna { input:
            vcf = gatk.snp_vcf,
            bam = merged.bam,
            num_cpus = run_3d_dna_num_cpus,
            ram_gb = run_3d_dna_ram_gb,
            disk_size_gb = run_3d_dna_disk_size_gb,
            runtime_environment = runtime_environment,
        }
    }
}

task create_fasta_index {
    input {
        File fasta
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        gzip -dc ~{fasta} > ~{basename(fasta, ".gz")}
        samtools faidx ~{basename(fasta, ".gz")}
    >>>

    output {
        File fasta_index = "~{basename(fasta, '.gz')}.fai"
    }

    runtime {
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task create_gatk_references {
    input {
        File reference_fasta
        File reference_fasta_index
        String output_stem
        RuntimeEnvironment runtime_environment
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
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
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
        Int ram_gb = 128
        Int disk_size_gb = 1000
        RuntimeEnvironment runtime_environment
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
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}


task run_3d_dna {
    input {
        File vcf
        File bam
        Int num_cpus = 8
        Int disk_size_gb = 2000
        Int ram_gb = 100
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        export VCF_FILENAME=~{basename(vcf, ".gz")}
        gzip -dc ~{vcf} > ${VCF_FILENAME}
        bash \
            /opt/3d-dna/phase/run-hic-phaser-encode.sh \
            --threads ~{num_cpus} \
            --to-stage update_vcf \
            ${VCF_FILENAME} \
            ~{bam}
        gzip -n *.txt *.vcf *.assembly
        ls
    >>>

    output {
        File snp_vcf = "snp.out.vcf.gz"
        File hic_vcf = "snp.out_HiC.vcf.gz"

        # .hic files
        File hic_in= "snp.out.in.hic"
        File hic = "snp.out.out.hic"

        # Scaffold boundary files (Juicebox 2D annotation format)
        File scaffold_track = "snp.out.out_asm.scaffold_track.txt.gz"
        File superscaffold_track = "snp.out.out_asm.superscaf_track.txt.gz"
        File scaffold_track_in = "snp.out.in_asm.scaffold_track.txt.gz"
        File superscaffold_track_in = "snp.out.in_asm.superscaf_track.txt.gz"

        File assembly_in = "snp.out.in.assembly.gz"
        File assembly = "snp.out.out.assembly.gz"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory: "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
