version 1.0

struct RuntimeEnvironment {
    String docker
    String singularity
}

workflow make_restriction_site_locations {
    meta {
        version: "1.14.2"
        caper_docker: "encodedcc/hic-pipeline:1.14.2"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.14.2"
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
        String docker = "encodedcc/hic-pipeline:1.14.2"
        String singularity = "docker://encodedcc/hic-pipeline:1.14.2"
    }


    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }

    call make_restriction_site_locations_ { input:
            reference_fasta = reference_fasta,
            assembly_name = assembly_name,
            restriction_enzyme = restriction_enzyme,
            runtime_environment = runtime_environment,
        }
    }

task make_restriction_site_locations_ {
    input {
        File reference_fasta
        String assembly_name
        String restriction_enzyme
        RuntimeEnvironment runtime_environment
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
        disks: "local-disk 10 HDD"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
