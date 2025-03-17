# Configuration options

A description of every valid option in `config.yaml`.


```
individuals:
  ind_1:
    deepvariant_vcf: /path/to/vcf/deepvariant_phased_ind_1.vcf.gz
    condition1:
      label_id1:
      - /path/to/flnc/invid_1/condition1/IsoSeqX.flnc.bam
    condition2:
      label_id2:
        - /path/to/flnc/invid_1/condition2/IsoSeqX.flnc.bam
  ind_2:
    deepvariant_vcf: /path/to/vcf/indiv_2.haplotagged.vcf.gz
    condition1:
      label_id1:
        - /path/to/flnc/ind_2/condition1/IsoSeqX.flnc1.bam
        - /path/to/flnc/ind_2/condition1/IsoSeqX.flnc2.bam
        - /path/to/flnc/ind_2/condition1/IsoSeqX.flnc3.bam
    condition2:
      label_id2:
        - /path/to/flnc/ind_2/condition2/IsoSeqX.flnc1.bam
        - /path/to/flnc/ind_2/condition2/IsoSeqX.flnc2.bam
        - /path/to/flnc/ind_2/condition2/IsoSeqX.flnc3.bam
  ind_3:
    deepvariant_vcf: 
    condition1:
      label_id1:
        - /path/to/flnc/ind_3/condition1/IsoSeqX.flnc1.bam
        - /path/to/flnc/ind_3/condition1/IsoSeqX.flnc2.bam
    condition2:
      label_id2:
        - /path/to/flnc/ind_3/condition2/IsoSeqX.flnc1.bam
        - /path/to/flnc/ind_3/condition2/IsoSeqX.flnc2.bam
reference_genome: /path/to/fasta/genome/hg38.fa
pigeon_annot: /path/to/reference_annotfile/gencode.v46.annotation.gtf
threads: 8
max_reads: 7000000
```
