import "../../hic.wdl" as hic

workflow test_create_hic {
	File pairs_file
    File chrsz_
     
	call hic.create_hic as test_create_hic_task { input:
		pairs_file = pairs_file,
    	chrsz_ = chrsz_
    }
}