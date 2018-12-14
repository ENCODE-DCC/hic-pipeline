##Encode DCC Hi-C pipeline tester
##Author: Ana Cismaru(anacismaru@gmail.com)
import "../../hic.wdl" as hic

workflow test_dedup {
    File merged_sort
   
    call hic.dedup as test_dedup_task { input:
    	merged_sort = merged_sort
    }
}