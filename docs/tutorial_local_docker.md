# Tutorial for general UNIX computers with docker

1. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    $ cd hic-pipeline
    ```

2. Install [Caper](https://github.com/ENCODE-DCC/caper), requires Python > 3.4.1
    ```bash
    $ pip install caper
    ```
    
3. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/template_one.json 
    $ PIPELINE_METADATA=metadata.json
    $ caper run workflow/main_workflow/hic.wdl --use-docker -i ${INPUT} -m ${PIPELINE_METADATA}
    ```

4. Execution of the pipeline on your local machine should take less than 10 minutes (delay could result from internet connection speed and the need to download the Docker image that is several 100s Mb large). You will be able to find all outputs on `cromwell-executions/hic/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.
Due to the minimal size of the input files, the last two steps of the pipeline (HiCCUPs and Arrowhead) will fail to detect any loops or domains. 

5. See full specification for [input JSON file](input.md).
