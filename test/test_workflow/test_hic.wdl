##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "hic.wdl" as hic
import "test/test_utils.wdl" as utils

workflow test_hic{
    #Used in pipeline
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    
    File restriction_sites
    File chrsz
    File reference_index

    #Reference items
    Array[Array[File]] ref_collisions 
    Array[Array[File]] ref_collisions_low
    Array[Array[File]] ref_unmapped 
    Array[Array[File]] ref_mapq0 
    Array[Array[File]] ref_alignable 
    Array[Array[File]] ref_sort_file 
    
    Array[File]ref_merged_collisions
    Array[File] ref_merged_collisions_low 
    Array[File] ref_merged_unmapped
    Array[File] ref_merged_mapq0 
    Array[File] ref_merged_alignable
    
    #Merge sort refputs
    Array[File] ref_merged_sort_file 
    #Dedup refputs
    Array[File] ref_dedup 
    
    #Merge pairs file refputs
    File ref_merged_pairs
    #Create hic refputs
    File ref_hic
    #TADs output
    File ref_tads
    #HiCCUps output
    File ref_hiccups

    #QC refputs
    #File ref_align_qc
    
    call hic.hic as test{ input:
        fastq = fastq
    }
    
    #Align task outputs
    Array[Array[File]] collisions = test.out_collisions
    Array[Array[File]] collisions_low = test.out_collisions_low
    Array[Array[File]] unmapped = test.out_unmapped
    Array[Array[File]] mapq0 = test.out_mapq0
    Array[Array[File]] alignable = test.out_alignable
    Array[Array[File]] sort_file = test.out_sort_file
    
    #Merge task outputs
    Array[File] merged_collisions = test.out_merged_collisions
    Array[File] merged_collisions_low = test.out_merged_collisions_low
    Array[File] merged_unmapped = test.out_merged_unmapped
    Array[File] merged_mapq0 = test.out_merged_mapq0
    Array[File] merged_align = test.out_merged_align
    
    #Merge sort outputs
    Array[File] merge_sort = test.out_merge_sort
    #Dedup outputs
    Array[File] dedup = test.out_dedup
    
    #TODO FIGURE OUT HOW TO COMPARE THESE
    #Merge pairs file outputs
    #File merged_pairs = test.out_merged_pairs
    #Create hic outputs
    #File hic = test.out_hic
    #TADs output
    #File tads = test.out_tads
    #HiCCUps output
    #File hiccups = test.out_hiccups

    #QC outputs
    #Array[File] align_qc = test.out_align_qc

    #TODO: CHECK EVERYTHING DONE IN TEST TASKS
    # call utils.compare_md5sum as one_dimension { input:
    #     labels = ["HiC Files"]   #Array[String]
    #     files = flatten([merged_collisions,merged_collisions_low,merged_unmapped,merged_mapq0,merged_align,merge_sort,dedup])         #Array[File]
    #     ref_files = flatten([ref_merged_collisions,ref_merged_collisions_low,ref_merged_unmapped,ref_merged_mapq0,ref_merged_align,ref_merge_sort,ref_dedup])    #Array[File]
    # }
}