import "../../hic.wdl" as hic

workflow test_bam2pairs {
	File bam 		    # bam file to convert into pairs
    File chrsz          # chromosome sizes file

	call hic.bam2pairs as test { input:
	    bam_file = bam,
    	chrsz_ = chrsz
	}

	call tail_of_pairs as strip_header { input:
		pairs = test.out_file
	}

 	output {
        File pairs_no_header = strip_header.no_header
    } 
}


task tail_of_pairs{
    File pairs
    
    command{
        sed 1,5d ${pairs} > no_header.pairs
    } 

    output{
        File no_header = glob("no_header.pairs")[0]
    }
}