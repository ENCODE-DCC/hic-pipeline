# Input JSON

An input JSON file includes all input parameters and metadata for running pipelines. Items 1), 2) and 3) are mandatory. Items 4) and 5) are optional so that our pipeline will use default values if they are not defined. However, 

* Mandatory

1. Reference genome.
2. Input data file paths/URIs.
3. Adapters to be trimmed.

* Optional

4. Pipeline parameters.
5. Resource settings for jobs.



## Templates

We provide two template JSON files for both single ended and paired-end samples. We recommend to use one of these input JSON files instead of that used in the tutorial section. These template JSON files include all parameters of the pipeline with default values defined.

* [template](../examples/template_one.json) for a single library sample


Let us take a close look at the following template JSON. Comments are not allowed in a JSON file but we added some comments to help you understand each parameter.
```javascript
{
    ////////// 1) Input FASTQ files //////////
    "hic.fastq": [[[
        "test/test_data/merged_read1.fastq.gz",
        "test/test_data/merged_read2.fastq.gz"
    ]]],

    ////////// 2) Reference genome chromosome sizes//////////
    "hic.chrsz": "test/test_data/ce10_selected.chrom.sizes.tsv",
    
    ////////// 3) Restriction sites locations in the reference genome sequence //////////
    "hic.restriction_sites": "test/test_data/ce10_selected_MboI.txt",

    ////////// 4) Reference genome index //////////
    "thic.reference_index": "test/test_data/ce10_selected.tar.gz",
    
    ////////// 5) Ligation site sequence //////////
    "hic.ligation_site": "AGCTAGCT"

}
```

## Reference genome

In order to run the HiC pipeline you will need to specify the bwa index file prepared using a referemnce genome sequence. We recommend using reference files from the ENCODE portal to enasure comparability of the analysis results.

|reference file description|assembly|ENCODE portal link|
|-|-|-|
|bwa index|hg19|[link](https://www.encodeproject.org/files/ENCFF807MUK/)|
|genome fasta|hg19|[link](https://www.encodeproject.org/files/male.hg19/)|
|chromosome sizes|hg19|[link](https://www.encodeproject.org/files/male.hg19.chrom.sizes/)|
|bwa index|GRCh38|[link](https://www.encodeproject.org/files/ENCFF643CGH/)|
|genome fasta|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)|
|chromosome sizes|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/)|

