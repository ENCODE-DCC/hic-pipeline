version 1.0

workflow make_restriction_site_locations {
    meta {
        version: "0.5.0"
        caper_docker: "encodedcc/hic-pipeline:0.5.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:0.5.0"
    }

    parameter_meta {
        assembly_name: "Name of genome assembly"
        reference_fasta: "FASTA file for the genome of interest for which to generate sites file"
        restriction_enzyme: "Restriction enzyme to generate sites file with"
    }

    input {
        File reference_fasta
        String assembly_name
        String restriction_enzyme
    }

    call make_restriction_site_locations_ { input:
            reference_fasta = reference_fasta,
            assembly_name = assembly_name,
            restriction_enzyme = restriction_enzyme,
        }
    }

task make_restriction_site_locations_ {
    input {
        File reference_fasta
        String assembly_name
        String restriction_enzyme
    }

    command <<<
        python3 "$(which generate_site_positions.py)" ~{restriction_enzyme} ~{assembly_name} ~{reference_fasta}
        gzip -n "~{assembly_name}_~{restriction_enzyme}.txt"
    >>>

    output {
        File restriction_site_locations = "~{assembly_name}_~{restriction_enzyme}.txt.gz"
    }

    runtime {
        cpu : "1"
        memory: "500 MB"
        disks: "local-disk 10 SSD"
    }
}
