# Tutorial for general UNIX computers with singularity

1. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    $ cd hic-pipeline
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ wget https://github.com/broadinstitute/cromwell/releases/download/36/cromwell-36.jar
    $ chmod +rx cromwell-36.jar
    ```

3. CHECK YOUR SINGULARITY VERSION FIRST AND UPGRADE IT TO A VERSION `>=2.5.2` OR PIPELINE WILL NOT WORK CORRECTLY.
    ```bash
    $ singularity --version
    ```

4. Build the singularity image for the pipeline. The following pulls the pipeline docker image, and uses that to construct the singularity image. The image will be stored in `~/.singularity`.
    ```bash
    SINGULARITY_PULLFOLDER=~/.singularity singularity pull quay.io/encode-dcc/hic-pipeline:template
    ```

5. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/template_one.json
    $ PIPELINE_METADATA=metadata.json
    $ java -jar  -Dconfig.file=backends/backend.conf -Dbackend.default=singularity cromwell-36.jar run workflow/main_workflow/hic.wdl -o workflow_opts/singularity.json -i ${INPUT} -m ${PIPELINE_METADATA}
    ```

6. You will be able to find all outputs on `cromwell-executions/hic/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

7. See full specification for [input JSON file](input.md).
