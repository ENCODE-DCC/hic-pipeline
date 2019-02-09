# Tutorial for general UNIX computers with docker

1. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    $ cd hic-pipeline
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ wget https://github.com/broadinstitute/cromwell/releases/download/37/cromwell-37.jar
    $ chmod +rx cromwell-37.jar
    ```
    
3. Run a pipeline for the test sample. TO MODIFY
    ```bash
    $ INPUT=examples/template_one.json 
    $ PIPELINE_METADATA=metadata.json
    $ java -jar -Dconfig.file=backends/backend.conf cromwell-37.jar run workflow/main_workflow/hic.wdl -i ${INPUT} -o workflow_opts/docker.json -m ${PIPELINE_METADATA}
    ```

4. Execution of the pipeline on your local machine should take less than 10 minutes. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

5. See full specification for [input JSON file](input.md).
