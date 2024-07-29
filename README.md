# MAS-Seq processing pipeline

This is a pipeline designed to process MAS-Seq data from PacBio. The pipeline allows to combine different samples (PS0001, PS0002) within a condition (treated/untreated) per individual. 
Multiple individuals can be defined for a combined analysis. 

The pipeline takes flnc.bam files as input. We haplotag this flnc bam individually. Ideally, a phased vcf should be provided, but if no phased vcf is available, an artificial HP:i:0 tag will be assigned to the annotations. 

## Installation

You will need snakemake and the snakemake executor plugin to distribute jobs on a cluster. 

```
# Create the conda environment
conda create -n isoseq-smk python=3.11

# Activate the conda environment
conda activate isoseq-smk

# Install pip, snakemake, and snakemake-executor-plugin-slurm
conda install pip
pip install snakemake snakemake-executor-plugin-slurm
```

## Usage

First set up a configuration file. See `config/config.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.

We suggest using the following flags for execution control: 
```
-p (--printshellcmds): snakemake prints the shell commands that it executes for each rule.
-k (--keep-going): Snakemake continues executing the workflow even if some jobs fail.
--rerun-incomplete: Snakemake re-run incomplete
```
```
snakemake --profile profiles/slurm-executor/ --configfile config/config.yaml -p -k 
```


[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/Linting/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/black/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
