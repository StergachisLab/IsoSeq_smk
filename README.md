# MAS-Seq processing pipeline

This is a pipeline designed to process [MAS-Seq data from PacBio](https://isoseq.how/). In particular, this pipeline helps with the following steps: clustering, alignment, merging, collapsing and annotating with pigeon, as well as the inclusion of annotations back to the refined `flnc.bam` files. The pipeline allows to combine different samples (PS0001, PS0002) within a condition (treated/untreated) per individual.
Multiple individuals can be included for a combined analysis, ensuring that all samples receive consistent transcript/isoform IDs for downstream analyses.

The pipeline takes unaligned `flnc.bam` files as input. If multiple `flnc.bam` files exist for a single condition (replicates), they can all be added to the YAML file and will be merged prior to any processing. We highly recommend including a phased VCF file generated from HiPhase + DeepVariant (Fiber-seq) since we use this file to help haplotag each read. Ideally, this phased VCF should be specified on the `config.yaml` file. If no phased VCF is available, an artificial `HP:i:0` tag will be assigned to the annotations.

In this pipeline, we also incorporate custom Iso-Seq tags into the aligned `flnc.bam` files as part of the final output. These tags aim to integrate information derived from the Pigeon output file `pigeon_classification.txt`:

```
haplotype (HP:i:)
isoform (in:Z:)
structural_category (st:Z:)
associated_gene (gn:Z:)
associated_transcript (tn:Z:)
subcategory (sb:Z:)
isoform_haplotype_noncyclo_counts (hn:i:)
isoform_haplotype_cyclo_counts (hc:i:)
isoform_noncyclo_counts (sn:i:)
isoform_cyclo_counts (sc:i:)
```

## Installation

You will need snakemake and the snakemake executor plugin to distribute jobs on a cluster.

```
# 1) Create the conda environment by using the provided yml file
# conda config --set channel_priority strict
# conda env create --file snakemake_env.yml

# 2) Activate the conda environment
conda activate snakemake_env

```

## Testing

First set up a configuration file. See `config/test_config.yaml` for a commented example.
Then run snakemake with the following command pointing to your configuration file.

To verify that the pipeline functions as expected, we can utilize the provided test data and configuration file available in the repository.
When running the test with a single sample, it is necessary to exclude the final rule, isoranker_analysis, as this step requires multiple samples to execute correctly. To achieve this, we use the following command:

```
snakemake --configfile config/test_config.yaml \
  --profile profiles/slurm-executor/ \
  -p -k --verbose \
  --directory results_test \
  --omit-from isoranker_analysis
```

## Execution

We suggest using the following flags for execution control:

```
-p (--printshellcmds): snakemake prints the shell commands that it executes for each rule.
-k (--keep-going): Snakemake continues executing the workflow even if some jobs fail.
--rerun-incomplete: Enables Snakemake to re-execute any previously incomplete rules.
```

```
nohup snakemake -s workflow/Snakefile \
--profile profiles/slurm-executor/ \
--configfile config/config.yaml \
-p -k --verbose --retries=1 \
--directory results > nohup.capture
```

### Sharing conda environments on your HPC

For Stergachis lab members be sure to include the following in your .bashrc if you want to use a pre-computed conda env.

```
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
export APPTAINER_CACHEDIR=/mmfs1/gscratch/stergachislab/snakemake-conda-envs/apptainer-cache
```
