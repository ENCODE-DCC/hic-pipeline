import "hic.wdl" as hic

workflow test_hic{
    #Used in pipeline
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[Array[File]]] input_bams = [] #[lib_id] MAKE 3D like fastqs
    Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
    Array[File] input_merged_sort = []
    Array[File] input_dedup_pairs = []
    File? input_pairs
    File restriction_sites
    File chrsz
    File reference_index

    #Reference
    File hic_ref
    
    call hic.hic as test{ input:
    fastq = fastq,
    input_bams = input_bams,
    input_sort_files = input_sort_files,
    input_merged_sort = input_merged_sort,
    input_dedup_pairs = input_dedup_pairs,
    input_pairs = input_pairs,
    restriction_sites = restriction_sites,
    chrsz = chrsz,
    reference_index = reference_index
    }
    File result = test.out_file

    call md5sum { input:
    labels = ["HiC Files"]
    files = [result]
    ref_files = [hic_ref]
    }
}