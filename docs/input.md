# Input JSON

An input JSON file includes all input parameters and metadata for running pipelines. Items 1), 2), 3), 4) and 5) are mandatory. Item 6) is optional so that our pipeline will use default values if it is not defined.

* Mandatory

1. Input FASTQ file pairs.
2. Reference genome bwa index.
3. Reference genome chromosome sizes.
4. Restriction site locations in the reference genome sequence.
5. Name of the restriction enzyme

* Optional

6. Pipeline parameters.

## Templates

We provide three template JSON files for processing of a single library with one or more sequencing runs and for multiple libraries.
* [template](../examples/template_one.json) for a single sequencing run from a single library
* [template](../examples/template_two.json) for two sequencing runs from a single library
* [template](../examples/template_three.json) for two libraries, each having a single sequencing run

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
    "hic.reference_index": "test/test_data/ce10_selected.tar.gz",
    
    ////////// 5) Ligation site sequence //////////
    "hic.restriction_enzyme": "MboI"

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

You will also need a restriction map file appropriate for the restriction enzyme and assembly. MboI and DpnII share the same restriction map because they have the same recognition site.

|restriction enzymes|assembly|ENCODE portal link|
|-|-|-|
|DpnII, MboI|GRCh38|[link](https://www.encodeproject.org/files/ENCFF246KDZ/)|
|HindIII|GRCh38|[link](https://www.encodeproject.org/files/ENCFF509VQM/)|
|DpnII, MboI|hg19|[link](https://www.encodeproject.org/files/ENCFF955ICX/)|
|HindIII|hg19|[link](https://www.encodeproject.org/files/ENCFF997LWB/)|
