# Usage Guide

## Contents

* [Introduction](usage.md#introduction)
* [Installation](usage.md#installation)
* [Running Workflows](usage.md#running-workflows)
* [Inspecting Outputs](usage.md#inspecting-outputs)
* [Supported Platforms](usage.md#inspecting-outputs)
* [Using Singularity](usage.md#using-singularity)

## Introduction

Here are some notes on running and using this pipeline. Using [Caper](https://github.com/ENCODE-DCC/caper) is the canonical, supported and official way to use ENCODE Uniform Processing Pipelines. The example below uses the command `caper run`, which is the simplest way to run a single pipeline. For running multiple pipelines in a production setting using `caper server` is recommended. To find details on setting up the server, refer to the [Caper documentation](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md#usage).

## Installation

1. Git clone this pipeline.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    ```

2. Install [Caper](https://github.com/ENCODE-DCC/caper), requires `java` >= 1.8 and `python` >= 3.6, Caper >= 0.8.2.1 is required to run the pipeline. Caper is a Python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ pip install caper  # use pip3 if it doesn't work
    ```

3. Follow [Caper's README](https://github.com/ENCODE-DCC/caper) carefully to configure it for your platform (local, cloud, cluster, etc.)
> **IMPORTANT**: Configure your Caper configuration file `~/.caper/default.conf` correctly for your platform.

## Running Workflows

Make sure you have properly installed the pipeline as described in the [installation instructions](usage.md#installation). Make sure to run the following commands from the root of the repository (i.e. `cd hic-pipeline` if you have not done so already).

1. Prepare the input JSON file. This file contains the user-specified files and parameters to run the pipeline with. Different examples of input JSON files are available [here](./reference.md#inputs). Details about the different input parameters are available [here](./reference.md#input-desciptions). Copy and paste the entirety of the following command into your terminal (uses [heredoc](https://tldp.org/LDP/abs/html/here-docs.html) syntax) and press enter/return to create a file called `input.json` pointing to the test data in this repo as pipeline input:

```bash
cat << EOF > input.json
{
  "hic.assembly_name": "ce10",
  "hic.chrsz": "tests/data/ce10_selected.chrom.sizes.tsv",
  "hic.fastq": [
    [
      [
        "tests/data/merged_read1.fastq.gz",
        "tests/data/merged_read2.fastq.gz"
      ]
    ]
  ],
  "hic.reference_index": "tests/data/ce10_selected.tar.gz",
  "hic.restriction_enzymes": [
    "MboI"
  ],
  "hic.restriction_sites": "tests/data/ce10_selected_MboI.txt.gz"
}
EOF
```

2. Run the pipeline using Caper. The `-m` flag is used to give a memorable name to the metadata JSON file the pipeline will produce once it is finished describing the run. More details about the metadata JSON can be found in the [Cromwell documentation](https://cromwell.readthedocs.io/en/stable/api/RESTAPI/#workflowmetadataresponse)

```bash
  $ caper run hic-pipeline.wdl -i input.json -m hic_testrun_metadata.json
```

## Inspecting Outputs

Rather than needing to dig through the highly nested Cromwell output directories or complex JSON metadata, [Croo](https://github.com/ENCODE-DCC/croo) can be used to generate a more legible HTML table of paths to outputs. To invoke `croo`, run the following, passing a Cromwell metadata JSON file as input:

```bash
  $ croo "${PATH_TO_METADATA_JSON}"
```

## Supported Platforms

This pipeline can be run on a variety of platforms via Caper. For a list of supported platforms, see [Caper's list of built-in backends](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md#built-in-backends). These include local machines, Google Cloud Platform, Amazon Web Services, and a selection of HPC clusters, namely Slurm, PBS, and SGE. Furthermore, Caper provides the [ability to use a custom backend](https://github.com/ENCODE-DCC/caper#running-pipelines-on-a-custom-backend), which can be useful in getting it to work with your particular cluster or cluster configuration.

## Using Singularity

Caper comes with built-in support for using Singularity containers instead of Docker with `--singularity` option. This is useful in HPC environments where Docker usage is restricted. See [Caper documentation](https://github.com/ENCODE-DCC/caper/blob/master/DETAILS.md) for more information.
