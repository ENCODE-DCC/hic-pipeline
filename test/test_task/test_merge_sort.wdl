##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "../../hic.wdl" as hic

workflow test_merge_sort {
    Array[File] sort_file
   
    call hic.merge_sort as test_merge_sort_task { input:
     sort_files_ = sort_file
    }
}