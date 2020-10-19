# Reference

## Contents

- [Reference](#reference)
  - [Contents](#contents)
  - [Inputs](#inputs)
    - [Entrypoints](#entrypoints)
      - [From fastqs](#from-fastqs)
      - [From aligned BAM](#from-aligned-bam)
      - [From processed libraries](#from-processed-libraries)
      - [From merged libraries](#from-merged-libraries)
      - [From hic file for loop and TAD calls](#from-hic-file-for-loop-and-tad-calls)
      - [Generating restriction site files](#generating-restriction-site-files)
    - [Input descriptions](#input-descriptions)
    - [Reference files](#reference-files)
  - [Outputs](#outputs)
    - [Output descriptions](#output-descriptions)

## Inputs

This pipeline has several different supported modes of operation. As such there are various entrypoints into the pipeline, each with their own set of relevant inputs. The various entrypoints are described in detail [here](#entrypoints), and the individual parameters are described [here](#input-descriptions). We recommend first determining which entrypoint you need then cross-referencing the relevant input descriptions.

### Entrypoints

Under each individual entrypoint the inputs for that entrypoint are listed. To run the pipeline using that particular entrypoint you need only specify the required inputs.

#### From fastqs

Runs the pipeline from the very beginning starting from `fastq` files. `read_groups` is useful for adding extra metadata into the BAM headers.

*Required inputs*
* `fastq`
* `restriction_enzymes`
* `restriction_sites`
* `chrsz`
* `reference_index`

*Optional inputs*
* `read_groups`

#### From aligned BAM

Runs the pipeline from starting from `bam` files, skipping `bwa` alignment.

*Required inputs*
* `bams`
* `ligation_counts`

#### From processed libraries

Runs the pipeline from starting with files from libraries that have already been processed by the pipeline. While `alignment_stats` and `library_stats` are optional it is highly recommended to specify them to obtain QC values for the merged libraries.

*Required inputs*
* `input_dedup_pairs`

*Optional inputs*
* `alignment_stats`
* `library_stats`

#### From merged libraries

Runs the pipeline starting with single Juicer-format text file.

*Required inputs*
* `input_pairs`

#### From hic file for loop and TAD calls

Runs the pipeline starting with a `.hic` file for loop and TAD calling.

*Required inputs*
* `input_hic`

#### Generating restriction site files

Runs the restriction sites file generation step only.

*Required inputs*
* `reference_fasta`
* `restriction_site_locations_only`, should be set to `true`

### Input descriptions

* `fastq` is a twice nested array of input fastqs. The outermost level corresponds to the biological replicates. Each biological replicate then contains one or more technical replicates. Each technical replicate is an array containing two paired fastq files. In the example below there are two biological replicates, therefore there are two elements in the outermost array. The first biological replicate has two technical replicates, so it has two arrays of fastqs. The second biological replicate only has one technical replicate so it only has one array of fastqs. In all cases the fastq arrays have two elements corresponding to read one and read two of the paired-end sequencing.
```json
"hic.fastq": [
  [
    [
      "biorep1_techrep1_R1.fastq.gz",
      "biorep1_techrep1_R2.fastq.gz"
    ],
    [
      "biorep1_techrep2_R1.fastq.gz",
      "biorep1_techrep2_R2.fastq.gz"
    ]
  ],
  [
    [
      "biorep2_techrep1_R1.fastq.gz",
      "biorep2_techrep1_R2.fastq.gz"
    ]
  ]
]
```
* `read_groups` is an optional array of strings to be inserted into the BAM as the read group (@RG), passed via `samtools addreplacerg` `-r` option. One per SE read/read pair with nested array structure mirroring the `fastq` input
* `restriction_enzymes` is an array of names containing the restriction enzyme(s) used to generate the Hi-C libraries. Currently only `MboI`, `HindIII`, `DpnII`, and `none` are supported. `none` is useful for libraries like DNAse produced using a non-specific cutter.
* `restriction_sites` is a text file containing cut sites for the given restriction enzyme. For supported enzymes you can generate this using the [reference building entrypoint](#generating-restriction-site-files). Note that if you need to generate a sites file for a multiple digest or for an unsupported enzyme you will need to edit this script and run it yourself: https://github.com/aidenlab/juicer/blob/encode/misc/generate_site_positions.py
* `chrsz` is a chromosome sizes file for the desired assembly. It is a tab-separated text file whose rows take the form `[chromosome][TAB][size]`. You can find these on the ENCODE portal for some human and mouse assemblies, see [reference files](#reference-files)
* `reference_index` is a pre-generated BWA index for the desired assembly. Depending on your assembly you may also be able to find these on the ENCODE portal, see [reference files](#reference-files)
* `bams` is a nested array of aligned, unfiltered BAM files, organized by `[replicate][library]`
* `ligation_counts` is an array of text files containing ligation counts for the `fastq` pair, organized by `[biological replicate][techincal replicate]`. These should be calculated from `fastqs` using the Juicer `countligations` script: https://github.com/aidenlab/juicer/blob/encode/CPU/common/
countligations.sh
* `input_pairs` is a text file containing the paired fragments to use to generate `.hic` contact maps, a detailed description of the file format can be found [here](https://github.com/aidenlab/juicer/wiki/Pre#long-format)
* `input_hic` is an input `.hic` file which will be used to call loops and domains
* `input_dedup_pairs` is an a array consisting of text files in the [Juicer pre long format](https://github.com/aidenlab/juicer/wiki/Pre#long-format) of paired fragments, one per library.
* `alignment_stats` is an array consisting of text files of alignment stats, one per library. Use is recommended but not required when merging libraries in order to calculate quality metrics on the merged libraries.
* `library_stats` is an array consisting of text files of library stats, one per library. Use is recommended but not required when merging libraries in order to calculate quality metrics on the merged libraries.
* `reference_fasta` is FASTA file for the genome of interest to be used for generating restriction site locations. For the output locations file to have a descriptive filename it is also recommended to specify the `assembly_name`
* `no_bam2pairs` is a boolean which if `true` results in skipping generating `.pairs` files, defaults to `false`
* `no_call_loops` is a boolean which if `true` results in skipping calling loops, defaults to `false`. Since the loop calling requires GPUs it is recommended to set to `true` if you do not
* `no_call_tads` is a boolean which if `true` skips calling domains with arrowhead, defaults to `false`
* `include_mapq0_reads` is a boolean which if `true` results in chimeric reads (3+ alignments) with one MAPQ 0 read being classified as normal paired reads with the MAPQ 0 read discarded. If `false` such reads will be classified as low mapq collisions. Defaults to `false`
* `cpu` is number of threads to use for `bwa` alignment, it is recommended to leave at the default value.
* `assembly_name` is name of assembly to insert into hic file header, recommended to specify for reproducibility otherwise the resulting `.hic.` file may have variable data in the header (the matrix contents will still be the same).

### Reference files

In order to run the pipeline from the beginning you will need to specify the `bwa` index and chromosome sizes file. We recommend using reference files from the ENCODE portal to ensure comparability of the analysis results. Links to the reference `fasta` files are included in case you need to generate a custom restriction sites file.

|reference file description|assembly|ENCODE portal link|
|-|-|-|
|bwa index|GRCh38|[link](https://www.encodeproject.org/files/ENCFF643CGH/)|
|genome fasta|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)|
|chromosome sizes|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/)|
|bwa index|hg19|[link](https://www.encodeproject.org/files/ENCFF807MUK/)|
|genome fasta|hg19|[link](https://www.encodeproject.org/files/male.hg19/)|
|chromosome sizes|hg19|[link](https://www.encodeproject.org/files/male.hg19.chrom.sizes/)|

In most cases you will also need a restriction map file appropriate for the restriction enzyme and assembly. `MboI` and `DpnII` share the same restriction map because they have the same recognition site. If you don't see your enzyme here you can generate a custom sites file, see [generating restriction site files](#generating-restriction-site-files).

|restriction enzymes|assembly|ENCODE portal link|
|-|-|-|
|DpnII, MboI|GRCh38|[link](https://www.encodeproject.org/files/ENCFF246KDZ/)|
|HindIII|GRCh38|[link](https://www.encodeproject.org/files/ENCFF509VQM/)|
|DpnII, MboI|hg19|[link](https://www.encodeproject.org/files/ENCFF955ICX/)|
|HindIII|hg19|[link](https://www.encodeproject.org/files/ENCFF997LWB/)|

## Outputs

The exact outputs of the pipeline will depend on the combination of inputs. When the pipeline completes check Caper/Cromwell metadata JSON file, particularly the top-level key `outputs` for values that are not `null`. Descriptions of the individual outputs are [below](#output-descriptions).

A draft document describing the pipeline outputs and quality control (QC) values is [available on the ENCODE portal](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf).

### Output descriptions

* `restriction_site_locations` is the generated restriction sites file
* `alignable_bam` is an array of filtered BAM files, one per biological replicate
* `out_pairs` is an array of files in [`pairs` format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md), one per biological replicate
* `out_dedup` is an array of files in [Juicer long format](https://github.com/aidenlab/juicer/wiki/Pre#long-format), one per biological replicate
* `library_complexity_stats_json` is an array of library complexity QC statistics in JSON format, one per biological replicate.
* `stats` is an array of library QC statistics in JSON format, one per biological replicate. It includes statistics describing the quantity and nature of the Hi-C contacts.
* `alignment_stats_` is an array of arrays of alignment QC statistics in plain text, one per technical replicate.
* `bams_with_read_group` is an array of arrays of unfiltered BAM files, one per technical replicate, it is only applicable if `read_groups` was specified as a pipeline input.
* `merged_stats_json` is a JSON file containing alignment and library statistics for merged libraries
* `out_hic_1` is a [`.hic` file](https://github.com/aidenlab/juicer/wiki/Data#hic-files) containing the contact matrix filtered by MAPQ >= 1
* `out_hic_30` is a [`.hic` file](https://github.com/aidenlab/juicer/wiki/Data#hic-files) containing the contact matrix filtered by MAPQ >= 30
* `out_tads` contains `arrowhead` domain calls in the Juicer format described [here](https://github.com/aidenlab/juicer/wiki/Arrowhead#domain-list-content)
