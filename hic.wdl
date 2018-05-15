workflow hic {
    Array[Array[File]] fastq_files
    File restriction_sites
    File chrsz
    File reference_index


    Int fastqs_len = length(fastq_files)
    
    scatter(i in range(fastqs_len)){
        call align { input:
            restriction = restriction_sites,
            fastqs = fastq_files[i],
            chrsz = chrsz,
            idx_tar = reference_index
        }
    }

    
    Array[Array[File]] bam_files = [
        align.collisions,
        align.collisions_low_mapq,
        align.unampped,
        align.mapq0,
        align.alignable
    ]

    Int bams_len = length(bam_files)
    scatter(i in range(bams_len)){
        call merge { input:
            bams = bam_files[i]
        }
    }


    call merge_sort { input:
        sort_files = align.sort_file,
        
    }

    # we can collect the alignabel.bam using the array merge.out_file
    #call dedup { input:
    #    merged_sort = merge_sort.out_file
    #}

    #call create_hic { input:
    #    chrsz = chrsz,
    #    pairs_file = dedup.out_file
    #}

    #call call_tads { input:
    #    hic_file = create_hic.out_file
    #}
}


task align {
	File idx_tar 		# reference bwa index tar
	Array[File] fastqs 	# [read_end_id]
    File chrsz          # chromosome sizes file
    File restriction    # restriction enzyme sites in the reference genome

    command {       
        mkdir data && cd data && mkdir fastq && mkdir reference
        data_path=$(pwd)
        cd fastq
        ln -s ${fastqs[0]} $(pwd)/frag_R1.fastq.gz
        ln -s ${fastqs[1]} $(pwd)/frag_R2.fastq.gz
        cd ../reference && tar -xvf ${idx_tar}
        index_folder=$(ls)
        cd $index_folder
        reference_fasta=$(ls | head -1) 
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ../..
        bash /opt/scripts/juicer.sh -D /opt -d $data_path -S alignonly -z $reference_index_path -p ${chrsz} -y ${restriction} -s MboI
    }

    output {
        File collisions = glob("data/splits/*_collisions.bam")[0]
        File collisions_low_mapq = glob("data/splits/*_collisions_low_mapq.bam")[0]
        File unampped = glob("data/splits/*_unmapped.bam")[0]
        File mapq0 = glob("data/splits/*_mapq0.bam")[0]
        File alignable = glob("data/splits/*_alignable.bam")[0]
        File sort_file = glob("data/splits/*.sort.txt")[0]
       
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05022018"
        cpu : 32
        memory: "64G"
    }
}

task merge {
    Array[File] bams

    command {
        samtools merge merged.bam ${sep=' ' bams}   
    }

    output {
        File out_file = glob('merged.bam')[0]
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05022018"
    }
}

task merge_sort {
    Array[File] sort_files

    command {
        sort -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n --parallel=8 -S 10% ${sep=' ' sort_files}  > merged_sort.txt
    }

    output {
        File out_file = glob('merged_sort.txt')[0]
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05022018"

        #> 8 processors
        #> a lot of memory
    }
}

task dedup {
    File merged_sort

    command {
        touch dups.txt
        touch optdups.txt
        touch merged_nodups.txt
        awk -f /opt/scripts/common/dups.awk ${merged_sort}
    }

    output {
        File out_file = glob('merged_nodups.txt')[0]
    }

    runtime {
        docker : "quay.io/gabdank/juicer:encode05022018"
    }
}

task create_hic {
    File pairs_file
    File chrsz

    command {
        /opt/scripts/common/juicer_tools pre -s inter.txt -g inter_hists.m -q 1 ${pairs_file} inter.hic ${chrsz}
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

task call_loops {
    File hic_file

    command {
       /opt/scripts/common/juicer_tools hiccups ${hic_file} contact_loops
    }

    output {
        Array[File] out_file = glob('contact_loops/*.*')
    }

    runtime {
        docker : "quay.io/gabdank/nvidia-juicer:v1.2"
        gpuType: "nvidia-tesla-k80" gpuCount: 2 zones: ["us-west1-b"]

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
			# if text file, read in contents
			if f.endswith('.qc') or f.endswith('.txt') or \
				f.endswith('.log') or f.endswith('.out'):
				with open(f,'r') as fp:
					contents = fp.read()
				with open(ref_f,'r') as fp:
					ref_contents = fp.read()
			else:
				contents = ''
				ref_contents = ''
			matched = md5==ref_md5
			result['tasks'].append(OrderedDict([
				('label', label),
				('match', matched),
				('md5sum', md5),
				('ref_md5sum', ref_md5),
				('basename', os.path.basename(f)),
				('ref_basename', os.path.basename(ref_f)),
				('contents', contents),
				('ref_contents', ref_contents),
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