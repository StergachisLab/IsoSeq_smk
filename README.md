![Build Status](https://github.com/StergachisLab/IsoSeq_smk/actions/workflows/main.yml/badge.svg)
![Downloads](https://img.shields.io/github/downloads/StergachisLab/IsoSeq_smk/total)
![Last Commit](https://img.shields.io/github/last-commit/StergachisLab/IsoSeq_smk)
![License](https://img.shields.io/github/license/StergachisLab/IsoSeq_smk)
![Snakemake](https://img.shields.io/badge/snakemake-compatible-brightgreen.svg?logo=snakemake&logoColor=white)
![Conda](https://img.shields.io/badge/conda-envs-green?logo=anaconda)
![Reproducible](https://img.shields.io/badge/reproducible-yes-brightgreen.svg)
![language](https://img.shields.io/github/languages/top/StergachisLab/IsoSeq_smk)
![Lint & Format]((https://github.com/StergachisLab/IsoSeq_smk/actions/workflows/lint.yml/badge.svg)



# MAS-Seq processing pipeline

This is a pipeline designed to process [MAS-Seq data from PacBio](https://isoseq.how/). In particular, this pipeline helps with the following steps: clustering, alignment, merging, collapsing and annotating with pigeon, as well as the inclusion of annotations back to the refined `flnc.bam` files. The pipeline allows to combine different samples (PS0001, PS0002) within a condition (treated/untreated) per individual.
Multiple individuals can be included for a combined analysis, ensuring that all samples receive consistent transcript/isoform IDs for downstream analyses.

The pipeline takes unaligned `flnc.bam` files as input. If multiple `flnc.bam` files exist for a single condition (replicates), they can all be added to the YAML file and will be merged prior to any processing. We highly recommend including a phased VCF file generated from  DeepVariant (Fiber-seq from same sample) + HiPhase,  since we use this file to help haplotag each read. Ideally, this phased VCF should be specified on the `config.yaml` file. If no phased VCF is available, an artificial `HP:i:0` tag will be assigned to the annotations.

In this pipeline, we also incorporate custom Iso-Seq tags into the aligned `flnc.bam` files as part of the final output. These tags aim to integrate information derived from the Pigeon output file `pigeon_classification.txt`:

```
haplotype (HP:i:)
isoform (in:Z:)
structural_category (sc:Z:)
associated_gene (gn:Z:)
associated_transcript (tn:Z:)
subcategory (sb:Z:)
isoform_haplotype_noncyclo_counts (hn:i:)
isoform_haplotype_cyclo_counts (hc:i:)
isoform_noncyclo_counts (sn:i:)
isoform_cyclo_counts (sc:i:)
```

As an additional feature, we have integrated [Isoranker](https://github.com/yhhc2/IsoRanker) to perform statistical analyses on our sample set [last rule of the pipeline]. While this step is optional, it requires additional supporting documents to function correctly. We have provided the necessary files corresponding to our current `config.yaml` list of samples in the `docs/` folder.

The path to this folder is specified as an additional flag in the `config.yaml` file. If you are running a new analysis, please ensure that your `config.yaml` and `docs/` folder contain all the necessary files and file paths required for the pipeline to run successfully.

## Installation

You will need snakemake and the snakemake executor plugin to distribute jobs on a cluster.

```
# 1) Create the conda environment by using the provided yml file
# conda config --set channel_priority strict
# conda env create --file snakemake_env.yml

# 2) Activate the conda environment
conda activate snakemake_env
conda config --set solver classic
```

## Testing

First set up a configuration file. 
See  `config/github.config.yaml` for a template. 
See `config/test_config.yaml` for a commented example when running a single condition.  
See `config/config.yaml`  for a full example processing multiple samples with different SMRTcells and and conditions.
The pipeline works when full paths are used to define all files and folders. 
Then run snakemake with the following command pointing to your own customized configuration file.

To verify that the pipeline functions as expected, we can utilize the provided test data and configuration file available in the repository.
When running the test with a single sample, it is necessary to exclude the final rule, isoranker_analysis, as this step requires multiple samples to execute correctly. To achieve this, we use the following command:

```
snakemake --configfile config/github.config.yaml \
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
