import "../../workflow/sub_workflow/process_library.wdl" as hic

workflow test_dedup {
    File merged_sort
    File restriction_sites
    String restriction_enzyme
    File alignable_bam

    Map[String, String] restriction_enzyme_to_site = read_map("workflow/sub_workflow/restriction_enzyme_to_site.tsv")
    String ligation_site = restriction_enzyme_to_site[restriction_enzyme]

    call hic.dedup as test_dedup_task { input:
        merged_sort = merged_sort,
        restriction_sites = restriction_sites,
        ligation_site = ligation_site,
        alignable_bam = alignable_bam
    }

    call strip_headers { input:
        bam = test_dedup_task.deduped_bam
    }

    output{
        File deduped = test_dedup_task.out_file
        File deduped_no_header = strip_headers.no_header
        File library_complexity_json = test_dedup_task.library_complexity_json
        File stats_json = test_dedup_task.stats_json
    }
}

task strip_headers{
    File bam

    #it messes up with compare_md5.py since all the files with stripped header are having the same name
    command {
        FILE=$(basename "${bam}" ".bam")
        samtools view -h ${bam} | samtools view - > $FILE.no_header.sam
    }
    
    output{
        File no_header = glob("*.no_header.sam")[0]
    }
}