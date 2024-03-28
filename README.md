# Iso-seq processing pipeline

This is a pipeline designed to process Iso-Seq data from PacBio HiFi reads.

## Usage

First set up a configuration file. See `config/config.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
```
snakemake --profile profile/compute --configfile config/your_config.yaml
```

### Dependencies
Requirements for executing the pipeline are the following conda packages:
```
snakemake>=7.32.0
python<=3.11
```
**Note:**
Snakemake is currently not compatible with Python >=3.12 due to a change in f-strings

An example install could look like this:
```
conda create -n snakemake -c conda-forge -c bioconda 'snakemake>=7.32' 'python=3.11'
```


[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/Linting/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/black/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
