[![CircleCI](https://circleci.com/gh/ENCODE-DCC/hic-pipeline/tree/dev.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/hic-pipeline/tree/dev)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

HiC uniform processing pipeline
===================================================

# Directories
* `workflow_opts/` : Workflow option files (`.json`)
* `examples/` : input JSON examples
* `docker/` : Dockerfiles to build images for running pipelines, including GPU-enabled image for running `hiccups`
* `tests/` : test code, including WDL tests, see the [developer docs](docs/development.md) for more details on running them.

# Installation and tutorial

This pipeline supports many cloud platforms and cluster engines. It supports `docker` and `singularity` to resolve complicated software dependencies for the pipeline. A tutorial-based instruction for each platform will be helpful to understand how to run pipelines.

* Cloud platform (CLI)
  * [Google Cloud Platform](docs/tutorial_google.md)
* Stanford HPC server (CLI)
  * [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)  *soon to be added*
* Local Linux computers (CLI)
  * [Local system with `singularity`](docs/tutorial_local_singularity.md)
  * [Local system with `docker`](docs/tutorial_local_docker.md)
* Cluster engines (CLI)
  * [SLURM](docs/tutorial_slurm_singularity.md) *[partial documentation]*

## Input JSON file

[Input JSON file specification](docs/input.md) *[partial documentation]*

## Output directories

[Output directory specification](docs/output.md) *soon to be added*
