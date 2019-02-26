import "../../workflow/sub_workflow/process_library.wdl" as hic

workflow test_dedup {
    File merged_sort
    File restriction_sites
    String ligation_site
    File alignable_bam

    call hic.dedup as test_dedup_task { input:
        merged_sort = merged_sort,
        restriction_sites = restriction_sites,
        ligation_site = ligation_site,
        alignable_bam = alignable_bam
    }

    output{
        File deduped = test_dedup_task.out_file
        File deduped_bam = test_dedup_task.deduped_bam
    }
}