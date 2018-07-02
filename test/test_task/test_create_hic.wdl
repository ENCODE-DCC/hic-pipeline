import "hic.wdl" as hic

workflow test_create_hic {
     File pe_file
     File chrsz

     File ref_hic
     
   
     call hic.create_hic as test { input:
     pairs_file = pe_file,
     chrsz_ = chrsz
     }
        
	File result = test.out_file

     call hic.compare_md5sum { input :
		labels = [
			'hic'
		],
    
		files = [
			result
		],
        
		ref_files = [
			ref_hic
		],
	}

}