##Encode DCC Hi-C pipeline
##Author: Ana Cismaru(anacismaru@gmail.com)

import "hic_sub.wdl" as sub
workflow hic {
   #User inputs
    Array[Array[Array[File]]] fastq = [] #[lib_id][fastq_id][read_end_id]
    Array[Array[Array[File]]] input_bams = [] #[lib_id] MAKE 3D like fastqs
    Array[Array[File]] input_sort_files = [] #[lib_id] 2d Array [lib[sort1, sirt2]]
    Array[File] input_merged_sort = []
    Array[File] input_dedup_pairs = []
    File? input_pairs

    
    File restriction_sites
    File chrsz
    File reference_index

    #determine range of scatter
    Int lib_length = if (length(fastq) > 0) then length(fastq)
    else if  (length(input_bams) > 0) then length(input_bams)
    else if (length(input_sort_files) > 0) then length(input_sort_files)
    else length(input_merged_sort)


    #call sub workflow to support input from multiple different libraries
    scatter(i in range(lib_length)){
        File? sub_ms #to deal with multiple entries
        
        call sub.hic_sub{ input:
        sub_fastq = if (length(fastq) > 0) then fastq[i] else [],
        sub_input_bams = if (length(input_bams) > 0) then input_bams[i] else [],
        sub_input_sort_files = if (length(input_bams) > 0) then input_sort_files[i] else [],
        sub_input_merged_sort = if length(input_merged_sort)>0 then input_merged_sort[i] else sub_ms,

        sub_restriction_sites = restriction_sites,
        sub_chrsz = chrsz,
        sub_reference_index =reference_index
        }
    }
   
    call merge_pairs_file{ input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else hic_sub.out_dedup
    }

    #flatten input pairs(merged and deduped)

    call create_hic { input:
    #fix
        pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
        chrsz_ = chrsz 
        
    }

}

task merge_pairs_file{
    Array[File] not_merged_pe
    command{
         sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' not_merged_pe}  > merged_pairs.txt
    }
    output{
        File out_file = glob('merged_pairs.txt')[0]
    }
    runtime {
       docker : "quay.io/gabdank/juicer:encode05022018"

       #> 8 processors
       #> a lot of memory
   }
}

task create_hic {
   File pairs_file
   File chrsz_

   command {
       /opt/scripts/common/juicer_tools pre -s inter.txt -g inter_hists.m -q 1 ${pairs_file} inter.hic ${chrsz_}
   }

  output {
       # add inter_30 stuff
       File out_file = glob('inter.hic')[0]
   }

  runtime {
       docker : "quay.io/gabdank/juicer:encode05022018"
   }
}

