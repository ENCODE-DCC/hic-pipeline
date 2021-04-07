[![CircleCI](https://circleci.com/gh/ENCODE-DCC/hic-pipeline/tree/dev.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/hic-pipeline/tree/dev)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![MIT License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

# ENCODE Hi-C uniform processing pipeline

## Overview

The [ENCODE](https://www.encodeproject.org/) pipeline for processing Hi-C data based on [Juicer](https://github.com/aidenlab/juicer)

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

## Usage

To verify your installation, you can run the following pipeline with a test data set by invoking the following command from the root of the cloned repository.

> Note: this will incur some cost when running in cloud environments.

```bash
$ caper run hic.wdl -i tests/functional/json/test_hic.json --docker
```

For detailed usage, see [usage](docs/usage.md)

## Inputs

See [inputs](docs/reference.md#inputs)


## Outputs

See [outputs](docs/reference.md#outputs)

## Contributing

We welcome comments, questions, suggestions, bug reports, feature requests, and pull requests (PRs). Please use one of the existing Github issue templates if applicable. When contributing code, please follow the [Developer Guidelines](docs/CONTRIBUTING.md#developer-guidelines).
