# Reference

# Contents

- [Reference](#reference)
- [Contents](#contents)
- [Workflows](#workflows)
  - [Hi-C](#hi-c)
    - [Inputs](#inputs)
      - [Entrypoints](#entrypoints)
        - [From fastqs](#from-fastqs)
        - [From hic file for annotations](#from-hic-file-for-annotations)
      - [Input descriptions](#input-descriptions)
      - [Reference files](#reference-files)
    - [Outputs](#outputs)
  - [Genophase](#genophase)
    - [Inputs](#inputs-1)
      - [Reference Files](#reference-files-1)
    - [Outputs](#outputs-1)
  - [Diploidify](#diploidify)
    - [Inputs](#inputs-2)
      - [Reference Files](#reference-files-2)
    - [Outputs](#outputs-2)
  - [Megamap](#megamap)
    - [Inputs](#inputs-3)
      - [Reference Files](#reference-files-3)
    - [Outputs](#outputs-3)
  - [Generating restriction site files](#generating-restriction-site-files)
- [Troubleshooting](#troubleshooting)
  - [Failure during create_hic](#failure-during-create_hic)
  - [Failure during `hic.normalize_assembly_name`](#failure-during-hicnormalize_assembly_name)
  - [Generic out of memory (OOM) issue](#generic-out-of-memory-oom-issue)
  - [Generic out of disk issue](#generic-out-of-disk-issue)
  - [Failure during localizer, hiccups, or arrowhead](#failure-during-localizer-hiccups-or-arrowhead)

# Workflows

[hic.wdl](../hic.wdl) is for starting with raw data (fastq) and producing `.hic` contact matrices and annotations (loops, topologically associating domains (TADs), subcompartments, stripes, loop domains, and eigenvectors, and in the case of intact data, DHS tracks). It is based on [Juicer](https://github.com/aidenlab/juicer)

[genophase.wdl](../genophase,wdl) takes in BAM files and produces phased variant calls. It is based on [hic2gatk](https://github.com/aidenlab/hic2gatk) and the [3D-DNA pipeline](https://github.com/aidenlab/3d-dna/).

[megamap.wdl](../megamap,wdl) takes in BAM files. It is based on [Juicer](https://github.com/aidenlab/juicer)

[diploidify.wdl](../diploidify,wdl) takes in BAM files and produces phased variant calls. It is based on [Juicer](https://github.com/aidenlab/juicer) [hic2gatk](https://github.com/aidenlab/hic2gatk) and the [3D-DNA pipeline](https://github.com/aidenlab/3d-dna/).

## Hi-C

This is the workflow contained in [hic.wdl](../hic.wdl)

### Inputs

This pipeline has several different supported modes of operation. As such there are various entrypoints into the pipeline, each with their own set of relevant inputs. The various entrypoints are described in detail [here](#entrypoints), and the individual parameters are described [here](#input-descriptions). We recommend first determining which entrypoint you need then cross-referencing the relevant input descriptions.

You can create an input JSON for running the pipeline end-to-end, i.e. from fastqs to loop and domain calls, on a Hi-C experiment from the ENCODE portal using the provided [input JSON generation script](../scripts/make_input_json_from_portal.py). Before running it install the requirements with `pip install -r requirements-scripts.txt`. To invoke it, you must at the minimum provide the accession of the experiment on the portal. See the script's help text for documentation of usage and options (`python scripts/make_input_json_from_portal.py --help`). Currently the script only supports experiments with one of `HindIII`, `DpnII`, and `MboI` as the library fragmentation method.

#### Entrypoints

Under each individual entrypoint the inputs for that entrypoint are listed. To run the pipeline using that particular entrypoint you need only specify the required inputs.

##### From fastqs

Runs the pipeline from the very beginning starting from `fastq` files. `read_groups` is useful for adding extra metadata into the BAM headers.

*Required inputs*
* `fastq`
* `restriction_enzymes`
* `restriction_sites`
* `chrsz`
* `reference_index`

If you are using a standard reference assembly like `GRCh38` or `mm10`, you should also specify it for the `assembly_name` input.

##### From hic file for annotations

Runs the pipeline starting with a `.hic` file for producing annotations.

*Required inputs*
* `input_hic`

#### Input descriptions

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
* `assembly_name` is name of assembly, defaults to "unknown". If the assembly is supported by Juicer Tools `pre` then `.hic` file creation will use Juicer Tools' internal chrom sizes instead of the inputted `chrsz`, see [`Pre` documentation](https://github.com/aidenlab/juicer/wiki/Pre#usage) for list of supported values. The pipeline does some normalization of this value internally, for instance `GRCh38` will be converted into the Juicer Tools-supported `hg38`.

The following flags offer control over which outputs are produced:
* `no_pairs` is a boolean which if `true` results in skipping generating `.pairs` files, defaults to `false`
* `no_call_loops` is a boolean which if `true` results in skipping calling loops, defaults to `false`. Since the loop calling requires GPUs it is recommended to set to `true` if you do not running the pipeline on Google Cloud
* `no_call_tads` is a boolean which if `true` skips calling domains with arrowhead, defaults to `false`
* `no_slice` is a boolean which if `true` skips calling subcompartments with SLICE, defaults to `false`
* `no_delta` is a boolean which if `true` results in skipping calling loops, domains, loop domains, and stripes with DELTA, defaults to `false`. Since DELTA requires GPUs it is recommended to set to `true` if you do not running the pipeline on Google Cloud
* `no_eigenvectors` is a boolean which if `true` skips producing A/B compartment bigWigs, defaults to `false`

The following parameters control resource usages. We recommend running with the default values, and only changing these if the workflow fails. You can find the defaults values by looking at the workflow-level `input` or the task-specific `input` in `hic.wdl`.

```
hic.align_num_cpus
hic.align_ram_gb_in_situ
hic.align_ram_gb_intact
hic.align_disk_size_gb_in_situ
hic.align_disk_size_gb_intact
hic.chimeric_sam_nonspecific_disk_size_gb
hic.chimeric_sam_specific_disk_size_gb
hic.dedup_ram_gb_in_situ
hic.dedup_ram_gb_intact
hic.dedup_disk_size_gb_in_situ
hic.dedup_disk_size_gb_intact
hic.create_hic_num_cpus
hic.create_hic_ram_gb
hic.create_hic_juicer_tools_heap_size_gb
hic.create_hic_disk_size_gb
hic.add_norm_num_cpus
hic.add_norm_ram_gb
hic.add_norm_disk_size_gb
hic.create_accessibility_track_ram_gb
hic.create_accessibility_track_disk_size_gb
```

#### Reference files

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

### Outputs

The exact outputs of the pipeline will depend on the combination of inputs. When the pipeline completes check Caper/Cromwell metadata JSON file, particularly the top-level key `outputs` for values that are not `null`. Descriptions of the individual outputs are [below](#output-descriptions).

A draft document describing the pipeline outputs and quality control (QC) values is [available on the ENCODE portal](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf).

## Genophase

This is the workflow contained in [genophase.wdl](../genophase,wdl). It calls SNPs using a GATK wrapper provided by `hic2gatk` and then phases the variants using `3d-dna`.

### Inputs

Three inputs are required to run the pipeline, `genophase.bams`, `genophase.gatk_bundle_tar`, and `genophase.reference_fasta`

The BAM files must have `RG` tags, which are required by GATK. A complete example input JSON is below.

```json
{
  "genophase.bams": [
    "https://www.encodeproject.org/files/ENCFF920JIL/@@download/ENCFF920JIL.bam",
    "https://www.encodeproject.org/files/ENCFF792QCT/@@download/ENCFF792QCT.bam"
  ],
  "genophase.gatk_bundle_tar": "https://www.encodeproject.org/files/ENCFF975DYE/@@download/ENCFF975DYE.tar.gz",
  "genophase.reference_fasta": "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
}
```

The following parameters can be used to adjust resources, however we recommend using the defaults on first pass. The default values are available in the WDL in each task input definition.

```
genophase.gatk_num_cpus
genophase.gatk_disk_size_gb
genophase.gatk_ram_gb
genophase.run_3d_dna_num_cpus
genophase.run_3d_dna_disk_size_gb
genophase.run_3d_dna_ram_gb
```

#### Reference Files

Pipeline uses reference variants for GATK. 4 VCFs and their indices are bundled into a tarball here: https://www.encodeproject.org/files/ENCFF975DYE/ . It contains Mills, Omni, Hapmap, and dbsnp variants from https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/

In addition, a reference fasta is required. It is downloadable from the portal [here](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)

### Outputs

The main output is the phased VCF, `run_3d_dna.hic_vcf` which has the filename "snp.out_HiC.vcf.gz".

Other useful outputs include the validation `.hic` file, `run_3d_dna.hic` which has the filename "snp.out.out.hic", and the `.psf` file `run_3d_dna.psf` which has the filename `out.psf`.

## Diploidify

This is the workflow contained in [diploidify.wdl](../diploidify,wdl). It takes the phased VCF file outputted by the `genophase` workflow, then applies it to BAM files from a run of the `hic.wdl` pipeline, resulting in the generation of haplotype-specific `.hic` contact maps and haplotype-specific DHS bigWig tracks.

### Inputs

The required inputs to the pipeline are `diploidify.bams`, `diploidify.chrom_sizes`, and `diploidify.vcf`.

`diploidify.vcf` is the `snp_hic` output of the `genophase` workflow. Do NOT use any other VCF from that pipeline, otherwise the results will be wrong.

```json
{
  "diploidify.bams": [
    "https://www.encodeproject.org/files/ENCFF920JIL/@@download/ENCFF920JIL.bam",
    "https://www.encodeproject.org/files/ENCFF792QCT/@@download/ENCFF792QCT.bam"
  ],
  "diploidify.chrom_sizes": "https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv",
  "diploidify.vcf": "gs://my/vcf/from/genophase/snp.out_HiC.vcf.gz"
}
```

Should you encounter any resource issues, the following parameters all control resource usage of various tasks:

```
diploidify.merge_num_cpus
diploidify.merge_ram_gb
diploidify.merge_disk_size_gb
diploidify.prepare_bam_num_cpus
diploidify.prepare_bam_ram_gb
diploidify.prepare_bam_disk_size_gb
diploidify.create_diploid_hic_num_cpus
diploidify.create_diploid_hic_ram_gb
diploidify.create_diploid_hic_disk_size_gb
diploidify.create_diploid_dhs_num_cpus
diploidify.create_diploid_dhs_ram_gb
diploidify.create_diploid_dhs_disk_size_gb
```

We recommend using the defaults first and only adjusting the resources if it fails.

#### Reference Files

A chromosome sizes file is required to run this pipeline. For GRCh38, the official ENCODE reference is [here](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv)

### Outputs

The outputs of this pipeline are all haplotype-specific. File from the two haplotypes are distinguished using `_r`/`_a` suffixes. This is just a convention and does not actually mean that one haplotype was used as the reference. It just means they are different haplotypes.

The two haplotype-specific `.hic` contact maps are available under the `hic_r` and `hic_a` keys in the output of the `create_diploid_hic` task.

The two haplotype-specific *raw* accessibility tracks are available under the `bigwig_raw_r` and `bigwig_raw_a` keys in the output of the `create_diploid_dhs` task.

The two haplotype-specific *corrected* accessibility tracks are available under the `bigwig_corrected_r` and `bigwig_corrected_a` keys in the output of the `create_diploid_dhs` task.

## Megamap

This is the workflow contained in [megamap.wdl](../megamap,wdl). It is very similar to the main `hic.wdl`. However, it starts from deduplicated bams from a previous run of the `hic.wdl` workflow. It also takes in the `.hic` files from the same previous runs as the bams in order to be able to calculate stats for the merged data.

### Inputs

The required inputs are `megamap.assembly_name`, `megamap.bams`, `megamap.chrom_sizes`, and `megamap.hic_files`. An example input is below:

```json
{
  "megamap.assembly_name": "GRCh38",
  "megamap.bams": [
    "https://www.encodeproject.org/files/ENCFF194KEQ/@@download/ENCFF194KEQ.bam",
    "https://www.encodeproject.org/files/ENCFF169MUQ/@@download/ENCFF169MUQ.bam"
  ],
  "megamap.chrom_sizes": "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
  "megamap.hic_files": [
    "https://www.encodeproject.org/files/ENCFF318JAP/@@download/ENCFF318JAP.hic",
    "https://www.encodeproject.org/files/ENCFF235LCO/@@download/ENCFF235LCO.hic"
  ]
}
```

#### Reference Files

A chromosome sizes file is required to run this pipeline. For GRCh38, the official ENCODE reference is [here](https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv)

### Outputs

The outputs are the similar to the main Hi-C pipeline. However, no bams are produced, just a merged `.hic` file and all the standard downstream annotations.

## Generating restriction site files

Use the WDL `make_restriction_site_locations.wdl` to generate the restriction sites file for your enzyme and reference fasta. Typically, you will want to just use the reference site locations file from the portal. However, if you must

*Required inputs*
* `reference_fasta`
* `restriction_enzyme`
* `assembly_name`

# Troubleshooting

If a workflow failed, you will see `Failed` status for that ID in `caper list`. To see what went wrong, run `caper debug WORKFLOW_ID`, replacing `WORKFLOW_ID` with the actual workflow ID given by `caper`. It will find the stdout/stderr of the failed tasks and print them to the terminal. Sometimes this information is not enough, and you may need to find the backend logs for the failed task, which capture issues that may have occurred before the task even started executing. The paths to these are available in the workflow metadata JSON which you can get via `caper metadata 1234`. Then, the paths to the logs are available under `.calls.[task_name].[shard_index].backendLogs`.

## Failure during create_hic

If you see `Killed` in the logs, try reducing the number of CPUs for `create_hic` task by setting `"hic.create_hic_num_cpus": 4` in your input JSON.

## Failure during `hic.normalize_assembly_name`

If this task fails, it almost certainly means something is wrong with your installation. Check that `docker` is installed and running if running locally.

## Generic out of memory (OOM) issue

If you `Killed` in the logs for the failed task, this means the machine ran out of RAM. Double the RAM for that failed task by updating the appropriate input (`hic.[task_name]_ram_gb`). Note that on Google Cloud the max allowed ram is ~624 GB, don't try more than that otherwise it will fail.

## Generic out of disk issue

If you see messages like `no space left on device`, this means that the task ran out of disk. Double the disk for that failed task by updating the appropriate input (`hic.[task_name]_disk_size_gb`).

## Failure during localizer, hiccups, or arrowhead

If you see things like "0 loops written to file", "data is too sparse", it means your data is not deep enough to call loops and/or domains. Turn off the loop and domain calls with `"hic.no_call_loops": true` and `"hic.no_call_tads": true` , then rerun the workflow.
