version 1.0

import "./hic.wdl"

workflow genophase {
    meta {
        version: "1.14.3"
        caper_docker: "encodedcc/hic-pipeline:1.14.3"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.14.3"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json"
    }

    input {
        File reference_fasta
        Array[File] bams
        # .tar.gz archive containing Ommi, Mills, Hapmap, and 1000G VCFs + indexes
        # https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
        # https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
        File? gatk_bundle_tar
        Int? gatk_num_cpus
        Int? gatk_disk_size_gb
        Int? gatk_ram_gb
        Int? run_3d_dna_num_cpus
        Int? run_3d_dna_disk_size_gb
        Int? run_3d_dna_ram_gb
        Boolean no_phasing = false

        String docker = "encodedcc/hic-pipeline:1.14.3"
        String singularity = "docker://encodedcc/hic-pipeline:1.14.3"
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
        bundle_tar = gatk_bundle_tar,
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

        call convert_psf_to_bedpe { input:
            psf = run_3d_dna.psf,
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
        File? bundle_tar
        Int num_cpus = 16
        Int ram_gb = 128
        Int disk_size_gb = 1000
        RuntimeEnvironment runtime_environment
    }

    String final_snp_vcf_name = "snp.out.vcf"
    String final_indel_vcf_name = "indel.out.vcf"

    command <<<
        mkdir bundle
        if [[ ~{if defined(bundle_tar) then "1" else "0"} -eq 1 ]]
        then
            tar -xvzf ~{bundle_tar} -C bundle
        fi
        mkdir reference
        mv ~{reference_fasta_index} ~{sequence_dictionary} ~{interval_list} reference
        gzip -dc ~{reference_fasta} > reference/~{basename(reference_fasta, ".gz")}
        run-gatk-after-juicer2.sh \
            -r reference/~{basename(reference_fasta, ".gz")} \
            ~{if defined(bundle_tar) then "--gatk-bundle bundle" else ""} \
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

        File psf = "out.psf"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory: "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task convert_psf_to_bedpe {
    input {
        File psf
        Int num_cpus = 1
        Int disk_size_gb = 1000
        Int ram_gb = 16
        RuntimeEnvironment runtime_environment
    }

    command <<<
        set -euo pipefail
        awk -f /opt/psf-to-bedpe/psf-to-bedpe.awk ~{psf} > "psf.bedpe"
    >>>

    output {
        File bedpe = "psf.bedpe"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk ~{disk_size_gb} HDD"
        memory: "~{ram_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
