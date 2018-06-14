import "hic_sub.wdl" as sub
workflow hic {
    #TODO: where do these actually belong? main or sub?
    File restriction_sites
    File chrsz
    File reference_index

    #Initialize variables to deal with multiple entry points
    # other input types (bam, nodup_bam)
    Array[File] input_bams = []  #[lib_id]
    Array[File] input_nodup_bams = [] #[lib_id]
    Array[File] input_sort_files = []

    #optional but important booleans
    Boolean align_only = false

    Array[Array[Array[File]]] fastq_libs = [] #[lib_id][fastq_id][read_end_id]
    Int lib_length = length(fastq_libs)
    
    #scatter within scatter to support different libraries
    scatter(i in range(lib_length)){
        call sub.hic_sub{ input:
        fastq_files = fastq_libs[i],
        sub_restriction_sites = restriction_sites,
        sub_chrsz = chrsz,
        sub_reference_index =reference_index

        }
    }
   
    call merge_pairs_file{ input:
        not_merged_pe = hic_sub.out_dedup
    }

    call create_hic { input:
        pairs_file = merge_pairs_file.out_file,
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

