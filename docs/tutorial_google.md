# Tutorial for Google Cloud Platform

All testing data is included in the repository. You won't have to download any data for testing our pipeline on Google Cloud.

1. Sign up for a Google account.
2. Go to [Google Project](https://console.developers.google.com/project) page and click "SIGN UP FOR FREE TRIAL" on the top left and agree to terms.
3. Set up a payment method and click "START MY FREE TRIAL".
4. Create a [Google Project](https://console.developers.google.com/project) `[YOUR_PROJECT_NAME]` and choose it on the top of the page.
5. Create a [Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser) `gs://[YOUR_BUCKET_NAME]` by clicking on a button "CREATE BUCKET" and create it to store pipeline outputs.
6. Find and enable following APIs in your [API Manager](https://console.developers.google.com/apis/library). Click a back button on your web brower after enabling each.
    * Compute Engine API
    * Google Cloud Storage (DO NOT click on "Create credentials")
    * Google Cloud Storage JSON API
    * Genomics API

7. Install [Google Cloud Platform SDK](https://cloud.google.com/sdk/downloads) and authenticate through it. You will be asked to enter verification keys. Get keys from the URLs they provide.
    ```bash
    $ gcloud auth login --no-launch-browser
    $ gcloud auth application-default login --no-launch-browser
    ```

8. If you see permission errors at runtime, then unset environment variable `GOOGLE_APPLICATION_CREDENTIALS` or add it to your BASH startup scripts (`$HOME/.bashrc` or `$HOME/.bash_profile`).
    ```bash
      unset GOOGLE_APPLICATION_CREDENTIALS
    ```

9. Set your default Google Cloud Project. Pipeline will provision instances on this project.
    ```bash
    $ gcloud config set project [YOUR_PROJECT_NAME]
    ```

10. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/hic-pipeline
    $ cd hic-pipeline
    ```

11. Download [cromwell](https://github.com/broadinstitute/cromwell). The pipeline has been tested with cromwell version 37.
    ```bash
    $ wget https://github.com/broadinstitute/cromwell/releases/download/37/cromwell-37.jar
    $ chmod +rx cromwell-37.jar
    ```

12. Run a pipeline using test dataset.
    ```bash
    $ PROJECT=[YOUR_PROJECT_NAME]
    $ BUCKET=gs://hic-pipeline-test-execution-bucket
    $ INPUT=examples/google_template.json 
    $ PIPELINE_METADATA=metadata.json

    $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=google -Dbackend.providers.google.config.project=${PROJECT} -Dbackend.providers.google.config.root=${BUCKET} cromwell-37.jar run workflow/main_workflow/hic.wdl -i ${INPUT} -o workflow_opts/docker.json -m ${PIPELINE_METADATA}
    ```

13. You will be able to find all outputs on your Google Cloud bucket. `gs://[YOUR_BUCKET_NAME] cromwell-executions/hic/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

14. See full specification for [input JSON file](input.md).
