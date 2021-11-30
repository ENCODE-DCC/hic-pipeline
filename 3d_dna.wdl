version 1.0

workflow phasing {
    meta {
        version: "1.4.0"
        caper_docker: "encodedcc/hic-pipeline:1.4.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.4.0"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/hic-pipeline/dev/croo_out_def.json"
    }

    input {
        File reference_fasta
        File bam
        String docker = "encodedcc/hic-pipeline:1.4.0"
    }

    call run_3d_dna { input:
        fasta = reference_fasta,
        merged_nodups = bam,
        docker = docker,
    }

    output {
        File hic = run_3d_dna.hic
    }
}

task run_3d_dna {
    input {
        File fasta
        File merged_nodups
        String docker
    }

    command <<<
        export MERGED_NODUPS_FILENAME=merged_nodups.txt
        gzip -dc ~{merged_nodups} > $MERGED_NODUPS_FILENAME
        bash /opt/run-3ddna-pipeline.sh ~{fasta} $MERGED_NODUPS_FILENAME
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
        docker: "~{docker}"
        memory: "32 GB"
    }
}
