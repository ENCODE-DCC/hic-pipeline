# SLURM

Using a generic SLURM cluster should be quite similar to Stanford Sherlock, which is a SLURM machine of a specific kind). The main differences are that Singularity installation may be different, please check with your cluster provider for installation instructions. You also will need to edit `workflow_opts/slurm.json` to include your information and the directories that contain your input data.

## Prerequisites

* Singularity version has to be `>=2.5.2`

* Have Java installed (installation method will depend on cluster)

* Check that the file `workflow_opts/slurm.json` is configured correctly. The `slurm.json` template file looks like this:

```
{
    "default_runtime_attributes" : {
      "slurm_partition": "[YOUR_SLURM_PARTITION]",
      "slurm_account": "[YOUR_SLURM_ACCOUNT]",
      "singularity_container" : "~/.singularity/hic-pipeline.simg",
      "singularity_command_options" : "--bind [DATA_DIR1],[DATA_DIR2],..."
    }
}
```

* Ensure your input data is located on the cluster

* Create an input JSON file specifying the paths to your input data and save it with a descriptive name, using a `.json` extension, to a memorable location, as you will need the path to this file to actually run the pipeline. You can find detailed specifications for this input file in the [input documentation](./input.md).

## Steps

1. Setup your partition, account and data:

Set your partition/account in `workflow_opts/slurm.json`. If your SLURM cluster does not require either a user's partition or account information, then remove them from this file. Otherwise, `YOUR_SLURM_PARTITON` or `YOUR_SLURM_ACCOUNT` will be used internally for `srun ... --partition YOUR_SLURM_PARTITION` and `srun ... --account YOUR_SLURM_PARTITION`, respectively. Replace `[DATA_DIR1],[DATA_DIR2],...` in the `singularity_command_options` with a comma-delimited list of paths to the directories that contain your input data.

2. Build the singularity image:

```bash
  $ SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://quay.io/encode-dcc/hic-pipeline:template
```

3. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/hic-pipeline
  $ cd hic-pipeline
```

4. Download Cromwell 35 using curl:

```bash
  $ curl -OL https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar
```

Alternatively, you can also use wget:

```bash
  $ wget -N -c https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar
```

5. Run the pipeline, replacing `INPUT_JSON_PATH` with the path to the input JSON file you created:

```bash
  $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=slurm_singularity cromwell-35.jar run hic-pipeline.wdl -i [INPUT_JSON_PATH] -o workflow_opts/slurm.json
```

You can also add ` -m [METADATA_FILE_PATH]` option to the above command, replacing `METADATA_FILE_PATH` with the path you would like Cromwell to write pipeline execution metadata to in JSON format, e.g. `my/metadata/path/my_experiment_metadata.json`. This is helpful for locating files produced by the pipeline and debugging.
