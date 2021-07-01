version 1.0

import "../../../hic.wdl" as hic

workflow test_create_eigenvector {
    input {
        File hic_file
        File chrom_sizes
    }

    call hic.create_eigenvector { input:
        hic_file = hic_file,
        chrom_sizes = chrom_sizes,
    }
}
