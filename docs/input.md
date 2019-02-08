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

* [template](../examples/template_se.json) for single ended sample
* [template](../examples/template_pe.json) for paired-end sample

Let us take a close look at the following template JSON. Comments are not allowed in a JSON file but we added some comments to help you understand each parameter.
```javascript
{
    ////////// 1) Reference genome //////////
    // Stanford servers: [GENOME]=hg38,hg19,mm10,mm9
    //   Sherlock: /home/groups/cherry/encode/pipeline_genome_data/[GENOME]_sherlock.tsv
    //   SCG4: /reference/ENCODE/pipeline_genome_data/[GENOME]_scg.tsv

    // Cloud platforms (Google Cloud, DNAnexus): [GENOME]=hg38,hg19,mm10,mm9
    //   Google Cloud: gs://encode-pipeline-genome-data/[GENOME]_google.tsv
    //   DNAnexus: dx://project-BKpvFg00VBPV975PgJ6Q03v6:data/pipeline-genome-data/[GENOME]_dx.tsv
    //   DNAnexus(Azure): dx://project-F6K911Q9xyfgJ36JFzv03Z5J:data/pipeline-genome-data/[GENOME]_dx_azure.tsv

    // On other computers download or build reference genome database and pick a TSV from [DEST_DIR].
    //   Downloader: ./genome/download_genome_data.sh [GENOME] [DEST_DIR]
    //   Builder (Conda required): ./conda/build_genome_data.sh [GENOME] [DEST_DIR]

    "atac.genome_tsv" : "/path_to_genome_data/hg38/hg38.tsv",

    ////////// 2) Input data files paths/URIs //////////

    // Read endedness
    "atac.paired_end" : true
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
|chromosome sizes|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/
)|

