# `config.yaml` Documentation  

This file configures input files, sample metadata, and reference genome paths. Please use absolute paths in this file. 

## Configuration Options  

### 1. `individuals`  

Defines sample information, including phased VCF files and long-read sequencing data. Condition labels are pre-established.

#### Structure:  
```yaml
individuals:
  <individual_id>:
    deepvariant_vcf: <path_to_vcf_file>
    <condition_name>:
      <label_id>:
        - <path_to_flnc_bam_file>
```

### 2. `reference_genome`
Path to the reference genome FASTA file.
```yaml
reference_genome: <path_to_fasta_file>
```

### 3. `pigeon_annot`
Path to the GTF annotation file used for transcript annotation.
Some tools like pigeon classify require a pre-indexed annotation file.
To index, please try manually creating the .pgi index before running the pipeline. 
```
pigeon index test_data/gtf/gencode.v46.annotation.gtf
```

```yaml
pigeon_annot: <path_to_gencode_gtf_file>
```

### 4. `docs_dir`
Path to the supporting isoranker documents.
Template uses docs nested inside our main IsoSeq_smk pipeline, with the relative path just above workflow. 
Feel free to specify an absolute path where your files are stored. File extensions do matter, please ensure all files have the proper extensions prior to submission.
```yaml
docs_dir: "../docs"
```

### 5. `threads`
Number of CPU threads to use for processing.
```yaml
threads: <num_threads>
```
Example: 

```yaml
individuals:
  ind_1:
    deepvariant_vcf: /path/to/vcf/deepvariant_phased_ind_1.vcf.gz
    untreated:
      label_A:
        - /path/to/flnc/ind_1/condition1/IsoSeqX.flnc.bam
    treated:
      label_B:
        - /path/to/flnc/ind_1/condition2/IsoSeqX.flnc.bam
  ind_2:
    deepvariant_vcf: /path/to/vcf/indiv_2.haplotagged.vcf.gz
    untreated:
      label_A:
        - /path/to/flnc/ind_2/condition1/IsoSeqX.flnc1.bam
        - /path/to/flnc/ind_2/condition1/IsoSeqX.flnc2.bam
        - /path/to/flnc/ind_2/condition1/IsoSeqX.flnc3.bam
    treated:
      label_B:
        - /path/to/flnc/ind_2/condition2/IsoSeqX.flnc.bam
reference_genome: /path/to/fasta/genome/hg38.fa
pigeon_annot: /path/to/reference_annotfile/gencode.v46.annotation.gtf
docs_dir: "../docs"
threads: 8
```
