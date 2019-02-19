[![CircleCI](https://circleci.com/gh/ENCODE-DCC/hic-pipeline/tree/dev-1.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/hic-pipeline/tree/dev-1)

HiC uniform processing pipeline
===================================================

# Directories
* `backends/` : Backend configuration files (`.conf`)
* `workflow_opts/` : Workflow option files (`.json`)
* `examples/` : input JSON examples
* `docker_image/` : Dockerfile
* `nvidia_docker_image/` : Dockerfile for Nvidia Docker: to use pipeline with GPU
* `test/` : test scripts for developers


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
