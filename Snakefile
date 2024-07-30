import os
import yaml
import glob

configfile: 'config/config.yaml'

# list of all individuals from the configuration file
individuals = config["individuals"].keys()
# list of conditions for each individual
conditions = ["treated", "untreated"]
# dictionary containing labels for each combination of individual and condition based on the configuration file
labels = {individual: {condition: list(config["individuals"][individual][condition].keys())
                       for condition in ["treated", "untreated"]}
          for individual in individuals}

# Define valid combinations of individual, condition, and label
combinations = [
    {"individual": ind, "condition": cond, "label": lbl}
    for ind in individuals
    for cond in ["treated", "untreated"]
    for lbl in labels[ind][cond]
]

# Restrict the values that wildcards can take from the config file
wildcard_constraints:
    individual = "|".join(individuals),
    condition = "treated|untreated",
    label = "|".join(set(label for ind in individuals for cond in conditions for label in labels[ind][cond]))

# check if "deepvariant_vcf" exists for an individual
def has_vcf(wc):
    return "deepvariant_vcf" in config["individuals"][wc.individual]

# get modified phased vcf file
def get_mod_phased_vcf(wildcards):
    phased_vcf = f"mod_vcf/{wildcards.individual}_mod.vcf.gz"
    if not has_vcf(wildcards):
        return f"mod_vcf/fake_{wildcards.individual}_mod.vcf.gz"
    return phased_vcf

# get the merged input for an individual, condition and label
def get_merge_input(wc):
    if wc.label in config["individuals"][wc.individual][wc.condition]:
        return config["individuals"][wc.individual][wc.condition][wc.label]
    else:
        raise KeyError(f"Label {wc.label} not found for individual {wc.individual} and condition {wc.condition}.")

# get vcf path for an individual
def get_vcf_path(wildcards):
    vcf_path = config['individuals'][wildcards.individual].get('deepvariant_vcf', None)
    if vcf_path:
        return vcf_path
    else:
        return f"mod_vcf/fake_{wildcards.individual}_mod.vcf.gz"

# combine labels for an individual
def combine_labels(wildcards):
    individual = wildcards.individual
    rtn = []
    for condition in ["treated", "untreated"]:
        for label in config["individuals"][individual][condition]:
            file = f"merged/{individual}_{condition}_{label}_labeled.bam"
            rtn.append(file)
    return rtn

# get haplotag files
def get_haplotag_files(wildcards):
    files = []
    for individual in config['individuals']:
        for condition in ["treated", "untreated"]:
            for label in config['individuals'][individual][condition]:
                files.append(f"whatshap/{individual}_{condition}_{label}.haplotagged.txt")
    return files

# get whatshap output files
def whatshap_outs(wc):
    template_f_path = "whatshap/{individual}_{condition}_{label}.tagged.bam"
    rtn = []
    for individual in config['individuals']:
        for condition in ["treated", "untreated"]:
            for label in config['individuals'][individual][condition]:
                rtn.append(template_f_path.format(
                    individual=individual,
                    condition=condition,
                    label=label
                ))
    return rtn 

rule all:
    input:
        expand("whatshap/{individual}_{condition}_{label}.tagged.bam",
               zip,
               individual=[comb["individual"] for comb in combinations],
               condition=[comb["condition"] for comb in combinations],
               label=[comb["label"] for comb in combinations]),
        "tag/dictionary_all.txt"


rule merge_individual_condition:
    conda:
        "envs/env.yml"
    input:
        get_merge_input
    output:
        bam = "merged/{individual}_{condition}_{label}_merged.bam"
    threads: config.get("threads", 40)
    benchmark: "benchmark/merge_individual_condition/{individual}_{condition}_{label}.txt"
    shell:
        "samtools merge -@ {threads} {output.bam} {input}"

rule label_reads_with_condition:
    conda:
        "envs/env.yml"
    input:
        bam = "merged/{individual}_{condition}_{label}_merged.bam"
    output:
        bam = "merged/{individual}_{condition}_{label}_labeled.bam"
    threads: config.get("threads", 40)
    benchmark: "benchmark/label_reads_condition/{individual}_{condition}_{label}.txt"
    resources:
        mem_mb=8000
    shell:
        r"""
        samtools view -@ {threads} -h {input.bam} | \
        awk -v id="{wildcards.individual}_{wildcards.condition}_{wildcards.label}" 'BEGIN {{OFS="\t"}} !/^@/ {{$1=id"_"$1; print}} /^@/ {{print}}' | \
        samtools view -bS - > {output.bam}
        """

rule align_sample:
    conda:
        "envs/env.yml"
    input:
        bam = "merged/{individual}_{condition}_{label}_labeled.bam"
    output:
        aligned = "aligned/{individual}_{condition}_{label}_aligned.bam"
    threads: 8
    benchmark: "benchmark/align_samples/{individual}_{condition}_{label}.txt"
    shell:
        """
        pbmm2 align {config[reference_genome]} {input.bam} {output.aligned} --preset ISOSEQ --sort -j {threads} --sort-memory 4G --log-level INFO
        """

rule modify_rg:
    conda:
        "envs/env.yml"
    input:
        bam = "aligned/{individual}_{condition}_{label}_aligned.bam"
    output:
        bam = "aligned/{individual}_{condition}_{label}_modified_aligned.bam",
        bai = "aligned/{individual}_{condition}_{label}_modified_aligned.bam.bai"
    shell:
        """
        samtools addreplacerg -r '@RG\\tID:rg1\\tSM:UnnamedSample\\tLB:lib1\\tPL:PACBIO' -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

rule modify_vcf:
    conda:
        "envs/env.yml"
    input:
        vcf = get_vcf_path
    output:
        mod_vcf = "mod_vcf/{individual}_mod.vcf.gz",
        mod_vcf_tbi = "mod_vcf/{individual}_mod.vcf.gz.tbi"
    threads: config.get("threads", 40)
    shell:
        r"""
        if [[ "{input.vcf}" =~ \.gz$ ]]; then
            zcat {input.vcf}
        else
            cat {input.vcf}
        fi | \
        bcftools view | \
        awk 'BEGIN {{FS=OFS="\t"}} /^#CHROM/ {{$10="UnnamedSample"; print; next}} {{print}}' | \
        bgzip -@ {threads} -c > {output.mod_vcf} && \
        tabix -p vcf {output.mod_vcf}
        """

rule fake_vcf:
    conda:
        "envs/env.yml"
    output:
        temp("mod_vcf/fake_{individual}_mod.vcf.gz"),
    shell:
        """
        touch {output}
        """

rule whatshap_stats:
    conda:
        "envs/env.yml"
    input:
        phased_vcf = "mod_vcf/{individual}_mod.vcf.gz"
    output:
        tsv = "whatshap/{individual}.phased_stats.tsv",
        block_list = "whatshap/{individual}.phased.blocks.txt",
        gtf = "whatshap/{individual}.phased.blocks.gtf"
    threads: config.get("threads", 40)
    shell:
        """
        whatshap stats --threads {threads} --tsv {output.tsv} --block-list {output.block_list} --gtf {output.gtf} {input.phased_vcf}
        """

rule whatshap_haplotag:
    conda:
        "envs/env.yml"
    input:
        phased_vcf_gz = get_mod_phased_vcf,
        bam = "aligned/{individual}_{condition}_{label}_modified_aligned.bam",
        reference = config['reference_genome']
    output:
        haplotagged_bam = "whatshap/{individual}_{condition}_{label}.haplotagged.bam",
        list_txt = "whatshap/{individual}_{condition}_{label}.list.txt"
    params:
        has_vcf = lambda wc: "true" if has_vcf(wc) else "false"
    threads: config.get("threads", 40)
    benchmark: "benchmark/whatshap_haplotag/{individual}_{condition}_{label}.txt"
    shell:
        """
        if [ {params.has_vcf} == "true" ]; then
            whatshap haplotag {input.phased_vcf_gz} {input.bam} -o {output.haplotagged_bam} --reference {input.reference} --output-haplotag-list {output.list_txt}
        else
            # add back fake haplotpag script
            ./scripts/add_hp0.sh {input.bam} {output.haplotagged_bam}  
        fi
        """

rule extract_read_info:
    conda:
        "envs/env.yml"
    input:
        bam = "whatshap/{individual}_{condition}_{label}.haplotagged.bam"
    output:
        txt = "whatshap/{individual}_{condition}_{label}.haplotagged.txt"
    benchmark: "benchmark/extract_read_info/{individual}_{condition}_{label}.txt"
    shell:
        """
        python ./scripts/extract_read_info.py {input.bam} {output.txt}
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
    benchmark: "benchmark/merge_samples/{individual}_all_conditions.txt"
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
        clustered = "clustered/{individual}_clustered.bam"
    threads: config.get("threads", 40)
    benchmark: "benchmark/cluster/{individual}_all_conditions.txt"
    shell:
        "isoseq3 cluster {input.bam} {output.clustered} --verbose -j {threads} --singletons"

rule align:
    conda:
        "envs/env.yml"
    input:
        "clustered/{individual}_clustered.bam"
    output:
        aligned="aligned/{individual}_aligned.bam"
    threads: 8
    benchmark: "benchmark/align/{individual}_aligned.txt"
    shell:
        "pbmm2 align {config[reference_genome]} {input} {output.aligned} --preset ISOSEQ --sort -j {threads} --sort-memory 4G --log-level INFO"

rule label_transcripts:
    conda:
        "envs/env.yml"
    input:
        bam = "aligned/{individual}_aligned.bam"
    output:
        bam = "aligned/{individual}_aligned_labeled.bam"
    threads: config.get("threads", 40)
    shell:
        r"""
        (
        samtools view -@ {threads} -h {input.bam} | awk -v id={wildcards.individual} 'BEGIN {{OFS="\t"}} !/^@/ {{$1=id"_"$1; print}} /^@/ {{print}}') | \
        samtools view -@ {threads} -bS - > {output.bam}
        """

rule merge_all_aligned:
    conda:
        "envs/env.yml"
    input:
        lambda wildcards: expand("aligned/{individual}_aligned_labeled.bam", individual=config['individuals'].keys())
    output:
        "aligned/all_individuals_aligned_merged.bam"
    threads: config.get("threads", 40)
    benchmark: "benchmark/merge_all_aligned/all_aligned_merged.txt"
    shell:
        """
        samtools merge -@ {threads} {output} {input}
        samtools index {output}
        """

rule collapse:
    conda:
        "envs/env.yml"
    input:
        "aligned/all_individuals_aligned_merged.bam"
    output:
        collapsed_gff="collapse/collapsed.gff",
        collapsed_abundance="collapse/collapsed.abundance.txt",
        collapsed_readstat="collapse/collapsed.read_stat.txt",
        collapsed_count="collapse/collapsed.flnc_count.txt"
    threads: config.get("threads", 40)
    benchmark: "benchmark/collapse/collapse.txt"
    shell:
        "isoseq3 collapse --num-threads {threads} {input} {output.collapsed_gff} --verbose"

rule pigeon_process:
    conda:
        "envs/env.yml"
    input:
        collapsed_gff="collapse/collapsed.gff",
        collapsed_count="collapse/collapsed.flnc_count.txt"
    output:
        sorted_gff="pigeon/sorted.gff",
        classification="pigeon/pigeon_classification.txt",
        report="pigeon/saturation.txt"
    params:
        pigeon_annot=config['pigeon_annot']
    benchmark: "benchmark/pigeon_process/pigeon.txt"
    shell:
        """
        pigeon sort {input.collapsed_gff} -o {output.sorted_gff}
        pigeon classify {output.sorted_gff} {params.pigeon_annot} {config[reference_genome]} --flnc {input.collapsed_count} -d pigeon -o 'pigeon'
        pigeon filter {output.classification}
        pigeon report pigeon/pigeon_classification.filtered_lite_classification.txt {output.report}
        """

rule create_dictionary:
    conda:
        "envs/r-base.yml"
    input:
        collapsed="collapse/collapsed.read_stat.txt",
        haplotag=get_haplotag_files,
        classification="pigeon/pigeon_classification.txt"
    output:
        dictionary="tag/dictionary_all.txt"
    benchmark: "benchmark/create_dictionary/dictionary.txt"
    shell:
        """
        Rscript ./scripts/create_dictionary.R {input.collapsed} {input.haplotag} {input.classification}
        """

rule add_tags_to_bam:
    conda:
        "envs/env.yml"
    input:
        bam="whatshap/{individual}_{condition}_{label}.haplotagged.bam",
        dictionary="tag/dictionary_all.txt"
    output:
        bam="whatshap/{individual}_{condition}_{label}.tagged.bam",
        bai="whatshap/{individual}_{condition}_{label}.tagged.bam.bai"
    benchmark: "benchmark/add_tags_to_bam/{individual}_{condition}_{label}_tagged.txt"
    shell:
        """
        ./scripts/add_tags_to_bam.sh {input.bam} {input.dictionary} {output.bam}
        """
