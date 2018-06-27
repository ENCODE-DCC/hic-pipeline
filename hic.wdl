##Encode DCC Hi-C pipeline
##Author: Ana Cismaru(anacismaru@gmail.com)

import "hic_sub.wdl" as sub
workflow hic {
   #User inputs
   # fastq files should look like filename_R1.fastq and filename_R2.fastq 
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
    else if  (length(input_bams) > 0) then length(input_bams) ##technically the number should be same for bams and sort_files
    else if (length(input_sort_files) > 0) then length(input_sort_files)
    else length(input_merged_sort)


    #call sub workflow to support input from multiple different libraries
    scatter(i in range(lib_length)){
        File? sub_ms #to deal with multiple entries
        
        call sub.hic_sub{ input:
        sub_fastq = if (length(fastq) > 0) then fastq[i] else [],
        sub_input_bams = if (length(input_bams) > 0) then input_bams[i] else [],
        sub_input_sort_files = if (length(input_sort_files) > 0) then input_sort_files[i] else [],
        sub_input_merged_sort = if length(input_merged_sort)>0 then input_merged_sort[i] else sub_ms,

        sub_restriction_sites = restriction_sites,
        sub_chrsz = chrsz,
        sub_reference_index =reference_index
        }
    }
   
    call merge_pairs_file{ input:
        not_merged_pe = if length(input_dedup_pairs)>0 then input_dedup_pairs else hic_sub.out_dedup
    }


    call create_hic { input:
        pairs_file = if defined(input_pairs) then input_pairs else merge_pairs_file.out_file,
        chrsz_ = chrsz     
    }
    
    #  call qc_report{ input:
    #  #TODO: FIND A WAY TO MERGE NORM_RES
    #  norm_res = norm.out_file,
    #  merged_nodups = merge_pairs_file.out_file,
    #  site_file = restriction_sites
    #  }
    
    # # call call_tads { input:
    # #   hic_file = create_hic.out_file
    # # }

   

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

task call_tads {
   File hic_file

   command {
       /opt/scripts/common/juicer_tools arrowhead ${hic_file} contact_domains --ignore_sparsity
   }

  output {
       Array[File] out_file = glob('contact_domains/*.bedpe')
   }

  runtime {
       docker : "quay.io/gabdank/juicer:encode05022018"
   }
}


task compare_md5sum {
	Array[String] labels
	Array[File] files
	Array[File] ref_files

	command <<<
		python <<CODE	
		from collections import OrderedDict
		import os
		import json
		import hashlib

		def md5sum(filename, blocksize=65536):
		    hash = hashlib.md5()
		    with open(filename, 'rb') as f:
		        for block in iter(lambda: f.read(blocksize), b""):
		            hash.update(block)
		    return hash.hexdigest()

		with open('${write_lines(labels)}','r') as fp:
			labels = fp.read().splitlines()
		with open('${write_lines(files)}','r') as fp:
			files = fp.read().splitlines()
		with open('${write_lines(ref_files)}','r') as fp:
			ref_files = fp.read().splitlines()

		result = OrderedDict()
		match = OrderedDict()
		match_overall = True

		result['tasks'] = []
		result['failed_task_labels'] = []
		result['succeeded_task_labels'] = []
		for i, label in enumerate(labels):
			f = files[i]
			ref_f = ref_files[i]
			md5 = md5sum(f)
			ref_md5 = md5sum(ref_f)
		
			matched = md5==ref_md5
			result['tasks'].append(OrderedDict([
				('label', label),
				('match', matched),
				('md5sum', md5),
				('ref_md5sum', ref_md5),
				('basename', os.path.basename(f)),
				('ref_basename', os.path.basename(ref_f)),
				]))
			match[label] = matched
			match_overall &= matched
			if matched:
				result['succeeded_task_labels'].append(label)
			else:
				result['failed_task_labels'].append(label)		
		result['match_overall'] = match_overall

		with open('result.json','w') as fp:
			fp.write(json.dumps(result, indent=4))
		match_tmp = []
		for key in match:
			val = match[key]
			match_tmp.append('{}\t{}'.format(key, val))
		with open('match.tsv','w') as fp:
			fp.writelines('\n'.join(match_tmp))
		with open('match_overall.txt','w') as fp:
			fp.write(str(match_overall))
		CODE
	>>>
	output {
		Map[String,String] match = read_map('match.tsv') # key:label, val:match
		Boolean match_overall = read_boolean('match_overall.txt')
		File json = glob('result.json')[0] # details (json file)
		String json_str = read_string('result.json') # details (string)
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}

task strip_headers{
    File bam
    
    command{
       samtools view ${bam} > no_header.sam
    } 
    output{
        File no_header = glob("no_header.sam")[0]
    }
    runtime {
        docker : "quay.io/gabdank/juicer:encode05232018"
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"     		
	}
}

# task qc_report{
#     File norm_res
#     File merged_nodups
#     File site_file
#     command{
#        ## Set ligation junction based on restriction enzyme
#       ##DO WE NEED THIS AND HOW DO I GET IT FROM RESTRICTION FILE
    
#       if [ -z "$ligation" ]
#           then
#                case $site in
#                HindIII) ligation="AAGCTAGCTT";;
#                DpnII) ligation="GATCGATC";;
#                MboI) ligation="GATCGATC";;
#                NcoI) ligation="CCATGCATGG";;
#                none) ligation="XXXX";;
#                *)  ligation="XXXX"
#                    echo "$site not listed as recognized enzyme. Using $site_file as site file"
#                    echo "Ligation junction is undefined"
#                esac
#            fi
          
#        #is this echo correct? It wanted the last line of the header file,but we did not make this file, but this should be the last line
#        #what is headfile??
#        echo "$0 $@"| awk '{printf"%-1000s\n", $0}' > inter_30.txt;  #what does this ; do?
# #Got rid of an $outputdir in front of the first inter_30.txt
# #Is merged_nodups coming from dedups? ligation?


#         cat ${norm_res} | awk -f /opt/scripts/common/stats_sub.awk >> inter_30.txt
#         java -cp /opt/scripts/common/ LibraryComplexity inter_30.txt >> inter_30.txt
#         /opt/scripts/common/statistics.pl -s ${site_file} -l $ligation -o inter_30.txt -q 30 ${merged_nodups}
    
        
#     }
#     output{
#         File out_file = glob("inter_30.txt")[0]
#     }
#     runtime{
#         docker : "quay.io/gabdank/juicer:encode05232018"
#         cpu : 32
#         memory: "64G"
#     }
# }





