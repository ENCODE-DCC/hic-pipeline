workflow arrowhead_test {
    File input_arrowhead
    
    call arrowhead { input:
        hic_file = input_arrowhead
    }
    
    output {
        File arrowhead_out = arrowhead.out_file
    } 
}

task arrowhead {
    
    File hic_file

    command {
        /opt/scripts/common/juicer_tools arrowhead ${hic_file} contact_domains 
    }

    output {
        File out_file = glob('contact_domains/*.bedpe')[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/hic-pipeline:template"
    }
}
