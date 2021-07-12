version 1.0

workflow create_hic_only {
    input {
        File pre
        File pre_index
        File stats
        File stats_hists
        File chrsz
        String assembly_name
        Int quality
        Int num_cpus
    }

    call create_hic { input:
        pre = pre,
        pre_index = pre_index,
        stats = stats,
        stats_hists = stats_hists,
        chrsz = chrsz,
        quality = quality,
        assembly_name = assembly_name,
        num_cpus = num_cpus,
    }
}

task create_hic {
    input {
        File pre
        File pre_index
        File stats
        File stats_hists
        Array[String] normalization_methods = []
        Int quality
        String? assembly_name
        File? chrsz
        File? restriction_sites
        Int num_cpus = 2
    }

    command <<<
        set -euo pipefail
        PRE_FILE=pre.txt
        PRE_INDEX_FILE=pre_index.txt
        RESTRICTION_SITES_FILENAME=restriction_sites.txt
        gzip -dc ~{pre} > $PRE_FILE
        gzip -dc ~{pre_index} > $PRE_INDEX_FILE
        ~{if defined(restriction_sites) then "gzip -dc " + restriction_sites + " > $RESTRICTION_SITES_FILENAME" else ""}
        # If the assembly name is empty, then we write chrsz path into file as usual, otherwise, use the assembly name instead of the path
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx20g \
            -jar /opt/scripts/common/juicer_tools.jar \
            pre \
            -n \
            ~{if defined(restriction_sites) then "-f $RESTRICTION_SITES_FILENAME" else ""} \
            -s ~{stats} \
            -g ~{stats_hists} \
            ~{if defined(assembly_name) then "-y " + assembly_name else ""} \
            -r 2500000,1000000 \
            -i $PRE_INDEX_FILE \
            --block-capacity 1000000 \
            --threads ~{num_cpus} \
            $PRE_FILE \
            inter_~{quality}.hic \
            ~{if defined(chrsz) then chrsz else assembly_name}
        java \
            -Ddevelopment=false \
            -Djava.awt.headless=true \
            -Xmx20g \
            -jar /opt/scripts/common/juicer_tools.jar \
            addNorm \
            ~{if length(normalization_methods) > 0 then "-k" else ""} ~{sep="," normalization_methods} \
            --threads ~{num_cpus} \
            inter_~{quality}.hic
    >>>

    output {
        File output_hic = "inter_~{quality}.hic"
    }

    runtime {
        cpu : "~{num_cpus}"
        disks: "local-disk 2000 SSD"
        memory : "30 GB"
    }
}
