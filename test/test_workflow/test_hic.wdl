import "../../workflow/main_workflow/hic.wdl" as hic_workflow

workflow test_hic {
    Array[Array[Array[File]]] fastq
    String ligation_site
    File restriction_sites
    File chrsz
    File reference_index
    File input_hic

    call hic_workflow.hic { input:
        fastq = fastq,
        ligation_site = ligation_site,
        restriction_sites = restriction_sites,
        chrsz = chrsz,
        reference_index = reference_index,
        input_hic = input_hic
    }

    call strip_headers { input:
        bam = hic.out_merged_align[0]
    }

    call strip_header { input:
        hic_file = hic.out_hic
    }

    output {
        File alignable_sam = strip_headers.no_header
        File alignable_pairs = hic.out_pairs[0]
        File merge_nodups = hic.out_dedup[0]
        File no_header_hic_file = strip_header.no_header

    }
}

task strip_headers{
    File bam

    #it messes up with compare_md5.py since all the files with stripped header are having the same name
    command {
        FILE=$(basename "${bam}" ".bam")
        samtools view -h ${bam} | samtools view - > $FILE.no_header.sam
    }
    
    output {
        File no_header = glob("*.no_header.sam")[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
	}
}

task strip_header {
    File hic_file

    command {
        hic_file=${hic_file}
        matrix_start=$(python3 /opt/straw/python/get_matrix_start.py $hic_file)
        tail -c +$matrix_start $hic_file > no_header.hic
        echo $(wc -c no_header.hic)
    }

    output {
        File no_header = glob("no_header.hic")[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
    }
}