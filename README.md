# Iso-seq processing pipeline

This is a pipeline designed to process Iso-Seq data from PacBio HiFi reads.

## Usage

First set up a configuration file. See `config/config.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
```
snakemake --configfile config/config.yaml -p --dag --notemp
```


[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/Linting/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/black/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
