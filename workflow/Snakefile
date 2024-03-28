# Import necessary libraries
import yaml
from os import path

# Load the configuration file
configfile: 'config.yaml'
wildcard_constraints: 
    individual = "|".join(config["individuals"]),
    conditions = "cyclo|non-cyclo",
    label = "PS\d+"

# get list of individuals on yml where vcf file is present
vcf_present =[]

for individual in config['individuals']:
    if not "deepvariant_vcf" in config['individuals'][individual]:
        continue
    file=config['individuals'][individual]["deepvariant_vcf"]
    if os.path.exists(file):
        vcf_present.append(individual)    

# Function to dynamically generate input paths for merging by condition
def get_sample_paths(individual, condition):
    print(individual, condition)
    x = config['individuals'][individual][condition]
    print(x)
    return config['individuals'][individual][condition]["path"]

# function to create sample label PS
def get_sample_labels_ind(individual, condition):
    return config['individuals'][individual][condition]["label"]

# function to combine samples from conditions cyclo/non-cyclo
def combine_labels(wildcards):
    individual = wildcards.individual
    rtn = []
    for condition in ["cyclo", "non-cyclo"]: 
        label = config["individuals"][individual][condition]["label"]
        file = f"merged/{label}_{individual}_{condition}_labeled.bam"
        rtn.append(file)
    return rtn

rule all:
    input:
        "pigeon/saturation.txt",
        expand("whatshap/{individual}.list.txt",individual=vcf_present)

# Rule to merge samples within each condition for an individual
rule merge_individual_condition:
    conda:
        "envs/env.yml"
    input:
        lambda wildcards: get_sample_paths(wildcards.individual, wildcards.condition)
    output:
        "merged/{label}_{individual}_{condition}_merged.bam"
    threads: config.get("threads", 40)
    shell:
        "samtools merge -@ {threads} {output} {input}"

rule label_reads_with_condition:
    conda:
        "envs/env.yml"
    input:
        "merged/{label}_{individual}_{condition}_merged.bam"
    output:
        bam = "merged/{label}_{individual}_{condition}_labeled.bam", 
        header = temp("{label}_{individual}_{condition}_header.sam")
    #params:
    #    label=lambda wildcards: wildcards.condition  # Directly use the condition name
    shell:
        """
        samtools view -H {input} > {output.header}
        samtools view {input} | \
        awk -v label="{wildcards.label}" 'BEGIN {{OFS="\\t"}} !/^@/ {{print label"_"$1, $0}} /^@/ {{print}}' | \
        cat {output.header} - | samtools view -bS > {output.bam}
        """

rule merge_samples:
    conda:
        "envs/env.yml"
    input:
        bam = combine_labels
    output:
        bam = "merged/{individual}_all_conditions_merged.bam",
        bai = "merged/{individual}_all_conditions_merged.bam.bai"
    threads: config.get("threads", 40)
    shell:
        """
        samtools merge -@ {threads} {output.bam} {input.bam}
        samtools index {output.bam}
        """

rule cluster:
    conda:
        "envs/env.yml"
    input:
        bam = "merged/{individual}_all_conditions_merged.bam",
        bai = "merged/{individual}_all_conditions_merged.bam.bai"
    output:
        clustered="clustered/{individual}_clustered.bam"
    threads: config.get("threads", 40)
    resources: 
        runtime = 60*24
    shell:
        "isoseq3 cluster {input.bam} {output.clustered} --verbose -j {threads} --singletons"

rule align:
    conda:
        "envs/env.yml"
    input:
        "clustered/{individual}_clustered.bam"
    output:
        aligned="aligned/{individual}_aligned.bam"
    threads: config.get("threads", 40)
    shell:
        "pbmm2 align {config[reference_genome]} {input} {output.aligned} --preset ISOSEQ --sort -j {threads} --sort-memory 40G --log-level INFO"

rule label_transcripts:
    conda:
        "envs/env.yml"
    input:
        "aligned/{individual}_aligned.bam"
    output:
        bam = "aligned/{individual}_aligned_labeled.bam",
        header = temp("{individual}_aligned_header.sam")
    #params:
    #    individual_label=lambda wildcards: config['individuals'][wildcards.individual]['label']
    shell:
        """
        samtools view -H {input} > {output.header}
        samtools view {input} | awk -v label={wildcards.individual} '{{print label"_"$$0}}' | cat {output.header} - | samtools view -bS > {output.bam}
        """

rule merge_all_aligned:
    conda:
        "envs/env.yml"
    input:
        lambda wildcards: expand("aligned/{individual}_aligned_labeled.bam", individual=config['individuals'].keys())
    output:
        "aligned/all_individuals_aligned_merged.bam"
    shell:
        """
        samtools merge -@ 40 {output} {input}
        samtools index {output}
        """

rule modify_vcf:
    conda:
        "envs/env.yml"
    input:
        vcf=lambda wildcards: config['individuals'][wildcards.individual]['deepvariant_vcf']
    output:
        mod_vcf="mod_vcf/{individual}_mod.vcf.gz",
        mod_vcf_tbi="mod_vcf/{individual}_mod.vcf.gz.tbi"
    shell:
        """
        if [[ "{input.vcf}" =~ \.gz$ ]]; then
            zcat {input.vcf}
        else
            cat {input.vcf}
        fi | \
        bcftools view | \
        awk 'BEGIN {{FS=OFS="\t"}} /^#CHROM/ {{$$10="UnnamedSample"; print; next}} {{print}}' | \
        bgzip -c > {output.mod_vcf} && \
        tabix -p vcf {output.mod_vcf}
        """

rule whatshap_stats:
    conda:
        "envs/env.yml"
    input:
        phased_vcf="mod_vcf/{individual}_mod.vcf.gz"
    output:
        tsv="whatshap/{individual}.phased_stats.tsv",
        block_list="whatshap/{individual}.phased.blocks.txt",
        gtf="whatshap/{individual}.phased.blocks.gtf"
    shell:
        "whatshap stats --tsv {output.tsv} --block-list {output.block_list} --gtf {output.gtf} {input.phased_vcf}"

rule whatshap_haplotag:
    conda:
        "envs/env.yml"
    input:
        phased_vcf_gz="mod_vcf/{individual}_mod.vcf.gz",
        bam="aligned/{individual}_aligned_labeled.bam",
        reference=config['reference_genome']
    output:
        haplotagged_bam="whatshap/{individual}.haplotagged.bam",
        list_txt="whatshap/{individual}.list.txt"
    shell:
        "whatshap haplotag {input.phased_vcf_gz} {input.bam} -o {output.haplotagged_bam} --reference {input.reference} --output-haplotag-list {output.list_txt}"

rule collapse:
    conda:
        "envs/env.yml"
    input:
        "aligned/all_individuals_aligned_merged.bam"
    output:
        collapsed_gff="collapse/collapsed.gff"
    shell:
        "isoseq3 collapse {input} {output.collapsed_gff} --verbose"

rule pigeon_process:
    conda:
        "envs/env.yml"
    input:
        collapsed_gff="collapse/collapsed.gff"
    output:
        report="pigeon/saturation.txt"
    params:
        pigeon_annot=config['pigeon_annot']
    shell:
        """
        pigeon sort {input.collapsed_gff} -o pigeon/sorted.gff
        pigeon classify pigeon/sorted.gff {params.pigeon_annot} {config[reference_genome]} --flnc collapse/collapsed.flnc_count.txt -d pigeon -o 'pigeon'
        pigeon filter pigeon/pigeon_classification.txt
        pigeon report pigeon/pigeon_classification.filtered_lite_classification.txt {output.report}
        """
