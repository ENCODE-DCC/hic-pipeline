version 1.0

import "./hic.wdl"

workflow diploidify {
    meta {
        version: "1.14.2"
        caper_docker: "encodedcc/hic-pipeline:1.14.2"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.14.2"
    }

    input {
        Array[File] bams
        # This is from genophase, snp.out_HiC.vcf.gz
        File vcf
        File chrom_sizes

        # Parameters
        Int quality = 30
        Array[Int] create_diploid_hic_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]

        # Resource params, specify to override the defaults
        Int? merge_num_cpus
        Int? merge_ram_gb
        Int? merge_disk_size_gb
        Int? prepare_bam_num_cpus
        Int? prepare_bam_ram_gb
        Int? prepare_bam_disk_size_gb
        Int? create_diploid_hic_num_cpus
        Int? create_diploid_hic_ram_gb
        Int? create_diploid_hic_disk_size_gb
        Int? create_diploid_dhs_num_cpus
        Int? create_diploid_dhs_ram_gb
        Int? create_diploid_dhs_disk_size_gb

        String docker = "encodedcc/hic-pipeline:1.14.2"
        String singularity = "docker://encodedcc/hic-pipeline:1.14.2"
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }

    call filter_chrom_sizes { input:
        chrom_sizes = chrom_sizes,
    }

    call hic.merge as merge { input:
        bams = bams,
        num_cpus = merge_num_cpus,
        ram_gb = merge_ram_gb,
        disk_size_gb = merge_disk_size_gb,
        runtime_environment = runtime_environment,
    }

    call prepare_bam { input:
        bam = merge.bam,
        quality = quality,
        vcf = vcf,
        chrom_sizes = chrom_sizes,
        num_cpus = prepare_bam_num_cpus,
        ram_gb = prepare_bam_ram_gb,
        disk_size_gb = prepare_bam_disk_size_gb,
    }

    call create_diploid_hic { input:
        bam = prepare_bam.filtered_bam,
        bam_index = prepare_bam.bam_index,
        vcf = vcf,
        chrom_sizes = filter_chrom_sizes.filtered_chrom_sizes,
        resolutions = create_diploid_hic_resolutions,
        num_cpus = create_diploid_hic_num_cpus,
        ram_gb = create_diploid_hic_ram_gb,
        disk_size_gb = create_diploid_hic_disk_size_gb,
    }

    call create_diploid_dhs { input:
        bam = prepare_bam.filtered_bam,
        bam_index = prepare_bam.bam_index,
        chrom_sizes = chrom_sizes,
        reads_to_homologs = create_diploid_hic.reads_to_homologs,
        num_cpus = create_diploid_dhs_num_cpus,
        ram_gb = create_diploid_dhs_ram_gb,
        disk_size_gb = create_diploid_dhs_disk_size_gb,
    }
}

task filter_chrom_sizes {
    input {
        File chrom_sizes
        String output_filename = "filtered.chrom.sizes"
    }

    command <<<
        python3 "$(which filter_chrom_sizes.py)" ~{chrom_sizes} ~{output_filename}
    >>>

    output {
        File filtered_chrom_sizes = output_filename
    }
}

task prepare_bam {
    input {
        File bam
        File vcf
        File chrom_sizes
        Int quality
        Int num_cpus = 8
        Int ram_gb = 64
        Int disk_size_gb = 2000
    }

    command <<<
        export CHROM_SIZES_FILENAME="assembly.chrom.sizes"
        mv ~{chrom_sizes} $CHROM_SIZES_FILENAME
        export VCF_FILENAME="snp.vcf"
        gzip -dc ~{vcf} > $VCF_FILENAME
        bash /opt/juicer/CPU/diploidify.sh \
            --from-stage prep \
            --to-stage prep \
            --vcf $VCF_FILENAME \
            --chrom-sizes $CHROM_SIZES_FILENAME \
            --mapq ~{quality} \
            --juicer-dir /opt \
            --phaser-dir /opt/3d-dna \
            ~{bam}
    >>>

    output {
        File filtered_bam = "reads.sorted.bam"
        File bam_index = "reads.sorted.bam.bai"
    }

    runtime {
        cpu: num_cpus
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}


task create_diploid_hic {
    input {
        File bam
        File bam_index
        File vcf
        File chrom_sizes
        Array[Int] resolutions
        Int num_cpus = 24
        Int ram_gb = 128
        Int disk_size_gb = 2000
    }

    command <<<
        mv ~{bam} "reads.sorted.bam"
        mv ~{bam_index} "reads.sorted.bam.bai"
        export CHROM_SIZES_FILENAME="assembly.chrom.sizes"
        mv ~{chrom_sizes} $CHROM_SIZES_FILENAME
        export VCF_FILENAME="snp.vcf"
        gzip -dc ~{vcf} > $VCF_FILENAME
        bash /opt/juicer/CPU/diploidify.sh \
            --from-stage hic \
            --to-stage hic \
            --vcf $VCF_FILENAME \
            --chrom-sizes $CHROM_SIZES_FILENAME \
            --resolutions ~{sep="," resolutions} \
            --threads-hic ~{num_cpus} \
            --juicer-dir /opt \
            --phaser-dir /opt/3d-dna
    >>>

    output {
        # r = reference, a = alternate haplotype but it's arbitrary
        File hic_r = "diploid_inter_r.hic"
        File hic_a = "diploid_inter_a.hic"
        File reads_to_homologs = "reads_to_homologs.txt"
    }

    runtime {
        cpu: num_cpus
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}

task create_diploid_dhs {
    input {
        File bam
        File bam_index
        File chrom_sizes
        File reads_to_homologs
        Int num_cpus = 2
        Int ram_gb = 128
        Int disk_size_gb = 1000
    }

    command <<<
        mv ~{bam} "reads.sorted.bam"
        mv ~{bam_index} "reads.sorted.bam.bai"
        export CHROM_SIZES_FILENAME="assembly.chrom.sizes"
        mv ~{chrom_sizes} $CHROM_SIZES_FILENAME
        bash /opt/juicer/CPU/diploidify.sh \
            --from-stage dhs \
            --to-stage dhs \
            --chrom-sizes $CHROM_SIZES_FILENAME \
            --reads-to-homologs ~{reads_to_homologs} \
            --juicer-dir /opt \
            --phaser-dir /opt/3d-dna
    >>>

    output {
        # r = reference, a = alternate haplotype but it's arbitrary
        File bigwig_raw_r = "diploid_inter_raw_r.bw"
        File bigwig_raw_a = "diploid_inter_raw_a.bw"
        File bigwig_corrected_r = "diploid_inter_corrected_r.bw"
        File bigwig_corrected_a = "diploid_inter_corrected_a.bw"
    }

    runtime {
        cpu: num_cpus
        memory: "~{ram_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}
