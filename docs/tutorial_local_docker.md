# Tutorial for general UNIX computers with docker

1. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    $ cd hic-pipeline
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell). The pipeline has been tested with cromwell version 37.
    ```bash
    $ wget https://github.com/broadinstitute/cromwell/releases/download/37/cromwell-37.jar
    $ chmod +rx cromwell-37.jar
    ```
    
3. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/template_one.json 
    $ PIPELINE_METADATA=metadata.json
    $ java -jar -Dconfig.file=backends/backend.conf cromwell-37.jar run workflow/main_workflow/hic.wdl -i ${INPUT} -o workflow_opts/docker.json -m ${PIPELINE_METADATA}
    ```

4. Execution of the pipeline on your local machine should take less than 10 minutes (delay could result from internet connection speed and the need to download the Docker image that is several 100s Mb large). You will be able to find all outputs on `cromwell-executions/hic/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.
Due to the minimal size of the input files, the last two steps of the pipeline (HiCCUPs and Arrowhead) will fail to detect any loops or domains. 

5. See full specification for [input JSON file](input.md).
