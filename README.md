# MAS-Seq processing pipeline

This is a pipeline designed to process MAS-Seq data from PacBio. The pipeline allows to combine different samples (PS0001, PS0002) within a condition (treated/untreated) per individual. 
Multiple individuals can be defined for a combined analysis. 

The pipeline takes flnc.bam files as input. We haplotag this flnc bam individually. Ideally, a phased vcf should be provided, but if no phased vcf is available, an artificial HP:i:0 tag will be assigned to the annotations. 


## Usage

First set up a configuration file. See `config/config.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
```
snakemake --configfile config/config.yaml -p
```


[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/Linting/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/black/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
