# Reference

## Contents

- [Reference](#reference)
  - [Contents](#contents)
  - [Workflows](#workflows)
  - [Inputs](#inputs)
    - [Entrypoints](#entrypoints)
      - [From fastqs](#from-fastqs)
      - [From hic file for annotations](#from-hic-file-for-annotations)
      - [Generating restriction site files](#generating-restriction-site-files)
    - [Input descriptions](#input-descriptions)
    - [Reference files](#reference-files)
  - [Outputs](#outputs)
  - [Genophase](#genophase)

## Workflows

[hic.wdl](../hic.wdl) is for starting with raw data (fastq) and producing `.hic` contact matrices and annotations (loops, TADs, subcompartments, stripes, loop domains, and A/B compartments). It is based on [Juicer](https://github.com/aidenlab/juicer)

[genophase.wdl](../genophase,wdl) takes in BAM files and produces phased variant calls. It is based on [hic2gatk](https://github.com/aidenlab/hic2gatk) and the [3D-DNA pipeline](https://github.com/aidenlab/3d-dna/).

`hic.wdl` is described below. For details on the other workflow, please jump to the [genophase](#genophase) section.

## Inputs

This pipeline has several different supported modes of operation. As such there are various entrypoints into the pipeline, each with their own set of relevant inputs. The various entrypoints are described in detail [here](#entrypoints), and the individual parameters are described [here](#input-descriptions). We recommend first determining which entrypoint you need then cross-referencing the relevant input descriptions.

You can create an input JSON for running the pipeline end-to-end, i.e. from fastqs to loop and domain calls, on a Hi-C experiment from the ENCODE portal using the provided [input JSON generation script](../scripts/make_input_json_from_portal.py). Before running it install the requirements with `pip install -r requirements-scripts.txt`. To invoke it, you must at the minimum provide the accession of the experiment on the portal. See the script's help text for documentation of usage and options (`python scripts/make_input_json_from_portal.py --help`). Currently the script only supports experiments with one of `HindIII`, `DpnII`, and `MboI` as the library fragmentation method.

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

If you are using a standard reference assembly like `GRCh38` or `mm10`, you should also specify it for the `assembly_name` input.

#### From hic file for annotations

Runs the pipeline starting with a `.hic` file for producing annotations.

*Required inputs*
* `input_hic`

#### Generating restriction site files

Use the WDL `make_restriction_site_locations.wdl` to generate the restriction sites file for your enzyme and reference fasta.

*Required inputs*
* `reference_fasta`
* `restriction_enzyme`
* `assembly_name`

### Input descriptions

* `fastq` is a twice nested array of input fastqs. The outermost level corresponds to the biological replicates. Each biological replicate then contains one or more `FastqPair` objects. In the example below there are two biological replicates, therefore there are two elements in the outermost array. The first biological replicate has two technical replicates, so it has two `FastqPair`s. The second biological replicate only has one technical replicate so it only has one `FastqPair`. The `FastqPair` can optionally contain a `read_group` which will be added to the aligned bam using `bwa mem` `-R` option. Note that this must contain the full header line, and that backslashes must be escaped in the JSON.
```json
"hic.fastq": [
  [
    {
      "read_1": "biorep1_techrep1_R1.fastq.gz",
      "read_2": "biorep1_techrep1_R2.fastq.gz"
    },
    {
      "read_1": "biorep1_techrep2_R1.fastq.gz",
      "read_2": "biorep1_techrep2_R2.fastq.gz"
    }
  ],
  [
    {
      "read_1": "biorep2_techrep1_R1.fastq.gz",
      "read_2": "biorep2_techrep1_R2.fastq.gz",
      "read_group": "@RG\\tID:myid"
    }
  ]
]
```
* `restriction_enzymes` is an array of names containing the restriction enzyme(s) used to generate the Hi-C libraries. Currently only `MboI`, `HindIII`, `DpnII`, and `none` are supported. `none` is useful for libraries like DNAse produced using a non-specific cutter.
* `ligation_site_regex` is a custom regular expression for counting ligation sites. If specified then `restriction_sites` file must be specified in the pipeline input. This can be just a single site, e.g. `ATGC`, or several sites wrapped in parentheses and separated by pipes, e.g. `(ATGC|CTAG)` (uses `grep -E` extended regular expression syntax)
* `restriction_sites` is a gzipped text file containing cut sites for the given restriction enzyme. For supported enzymes you can generate this using the [reference building entrypoint](#generating-restriction-site-files). Note that if you need to generate a sites file for a multiple digest or for an unsupported enzyme you will need to edit this script and run it yourself: https://github.com/aidenlab/juicer/blob/encode/misc/generate_site_positions.py
* `chrsz` is a chromosome sizes file for the desired assembly. It is a gzipped and tab-separated text file whose rows take the form `[chromosome][TAB][size]`. You can find these on the ENCODE portal for some human and mouse assemblies, see [reference files](#reference-files)
* `reference_index` is a pre-generated BWA index for the desired assembly. Depending on your assembly you may also be able to find these on the ENCODE portal, see [reference files](#reference-files)
* `input_hic` is an input `.hic` file which will be used to call loops and domains
* `normalization_methods` is an array of normalization methods to use for `.hic` file generation as per Juicer Tools `pre`. If not specified then will use `pre` defaults of `VC`, `VC_SQRT`, `KR`, and `SCALE`. Valid methods are `VC`, `VC_SQRT`, `KR`, `SCALE`, `GW_KR`, `GW_SCALE`, `GW_VC`, `INTER_KR`, `INTER_SCALE`, and `INTER_VC.
* `reference_fasta` is FASTA file for the genome of interest to be used for generating restriction site locations. For the output locations file to have a descriptive filename it is also recommended to specify the `assembly_name`
* `align_num_cpus` is number of threads to use for `bwa` alignment, it is recommended to leave at the default value.
* `create_hic_num_cpus` is number of threads to use for hic creation, it is recommended to leave at the default value. If you have an OOM error for Juicer Tools `pre`, which may occur for large experiments, then supply a small value such as `4`
* `assembly_name` is name of assembly, defaults to "unknown". If the assembly is supported by Juicer Tools `pre` then `.hic` file creation will use Juicer Tools' internal chrom sizes instead of the inputted `chrsz`, see [`Pre` documentation](https://github.com/aidenlab/juicer/wiki/Pre#usage) for list of supported values. The pipeline does some normalization of this value internally, for instance `GRCh38` will be converted into the Juicer Tools-supported `hg38`.

The following flags offer control over which outputs are produced:
* `no_pairs` is a boolean which if `true` results in skipping generating `.pairs` files, defaults to `false`
* `no_call_loops` is a boolean which if `true` results in skipping calling loops, defaults to `false`. Since the loop calling requires GPUs it is recommended to set to `true` if you do not running the pipeline on Google Cloud
* `no_call_tads` is a boolean which if `true` skips calling domains with arrowhead, defaults to `false`
* `no_slice` is a boolean which if `true` skips calling subcompartments with SLICE, defaults to `false`
* `no_delta` is a boolean which if `true` results in skipping calling loops, domains, loop domains, and stripes with DELTA, defaults to `false`. Since DELTA requires GPUs it is recommended to set to `true` if you do not running the pipeline on Google Cloud
* `no_eigenvectors` is a boolean which if `true` skips producing A/B compartment bigWigs, defaults to `false`

### Reference files

In order to run the pipeline from the beginning you will need to specify the `bwa` index and chromosome sizes file. We recommend using reference files from the ENCODE portal to ensure comparability of the analysis results. Links to the reference `fasta` files are included in case you need to generate a custom restriction sites file.

|reference file description|assembly|ENCODE portal link|
|-|-|-|
|bwa index|GRCh38|[link](https://www.encodeproject.org/files/ENCFF643CGH/)|
|genome fasta|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)|
|chromosome sizes|GRCh38|[link](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/)|
|bwa index|mm10|[link](https://www.encodeproject.org/files/ENCFF018NEO/)|
|chromosome sizes|mm10|[link](https://www.encodeproject.org/files/mm10_no_alt.chrom.sizes/)|

In most cases you will also need a restriction map file appropriate for the restriction enzyme and assembly. `MboI` and `DpnII` share the same restriction map because they have the same recognition site. If you don't see your enzyme here you can generate a custom sites file, see [generating restriction site files](#generating-restriction-site-files).

|restriction enzymes|assembly|ENCODE portal link|
|-|-|-|
|DpnII, MboI|GRCh38|[link](https://www.encodeproject.org/files/ENCFF246KDZ/)|
|HindIII|GRCh38|[link](https://www.encodeproject.org/files/ENCFF509VQM/)|
|DpnII, MboI|mm10|[link](https://www.encodeproject.org/files/ENCFF930KBK/)|

## Outputs

The exact outputs of the pipeline will depend on the combination of inputs. When the pipeline completes check Caper/Cromwell metadata JSON file, particularly the top-level key `outputs` for values that are not `null`. Descriptions of the individual outputs are [below](#output-descriptions).

A draft document describing the pipeline outputs and quality control (QC) values is [available on the ENCODE portal](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf).

## Genophase

