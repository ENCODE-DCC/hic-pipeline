# SLURM

Using a generic SLURM cluster should be quite similar to Stanford Sherlock, which is a SLURM machine of a specific kind). The main differences are that Singularity installation may be different, please check with your cluster provider for installation instructions. You also will need to edit `workflow_opts/slurm.json` to include your information and the directories that contain your input data.

## Prerequisites

* Singularity version has to be `>=2.5.2`

* Have Java and Python >= 3.6 installed (installation method will depend on cluster)

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

2. Get the code and move to the repo directory:

```bash
  $ git clone https://github.com/ENCODE-DCC/hic-pipeline
  $ cd hic-pipeline
```

3. Install [Caper](https://github.com/ENCODE-DCC/caper), requires Python >= 3.6. Caper >= 0.8.2.1 is required to run the pipeline.

```bash
  $ pip install caper
```

4. Run the pipeline, replacing `INPUT_JSON_PATH` with the path to the input JSON file you created:

```bash
  $ caper run hic-pipeline.wdl --use-singularity -b slurm -i [INPUT_JSON_PATH] -o workflow_opts/slurm.json
```

You can also add ` -m [METADATA_FILE_PATH]` option to the above command, replacing `METADATA_FILE_PATH` with the path you would like Caper to write pipeline execution metadata to in JSON format, e.g. `my/metadata/path/my_experiment_metadata.json`. This is helpful for locating files produced by the pipeline and debugging.
