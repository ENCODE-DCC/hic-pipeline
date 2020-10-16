# Tutorial for general UNIX computers with singularity

1. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    $ cd hic-pipeline
    ```

2. Install [Caper](https://github.com/ENCODE-DCC/caper), requires Python >= 3.6. Caper >= 0.8.2.1 is required to run the pipeline.
    ```bash
    $ pip install caper
    ```

3. CHECK YOUR SINGULARITY VERSION FIRST AND UPGRADE IT TO A VERSION `>=2.5.2` OR PIPELINE WILL NOT WORK CORRECTLY.
    ```bash
    $ singularity --version
    ```

4. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/template_one.json
    $ PIPELINE_METADATA=metadata.json
    $ caper run workflow/main_workflow/hic.wdl --use-singularity -i ${INPUT} -m ${PIPELINE_METADATA}
    ```

5. You will be able to find all outputs on `cromwell-executions/hic/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

6. See full specification for [input JSON file](input.md).
