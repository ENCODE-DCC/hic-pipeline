##Encode DCC Hi-C pipeline
import "test_hic.wdl" as hic

workflow test_multiple_libraries {
    Array[Array[Array[File]]] fastq
    String restriction_enzyme
    File restriction_sites
    File chrsz
    File reference_index

    Map[String, String] restriction_enzyme_to_site = read_map("workflow/sub_workflow/restriction_enzyme_to_site.tsv")
    String ligation_site = restriction_enzyme_to_site[restriction_enzyme]

    call hic.test_hic as test { input:
        fastq = fastq,
        chrsz = chrsz,
        restriction_sites = restriction_sites,
        reference_index = reference_index,
        restriction_enzyme = restriction_enzyme
    }

    output{
        File no_header_alignable_sam = test.no_header_alignable_sam
        File out_pairs = test.out_pairs
        File out_dedup = test.out_dedup
        File no_header_hic_1 = test.no_header_hic_1
        File no_header_hic_30 = test.no_header_hic_30
    }
}