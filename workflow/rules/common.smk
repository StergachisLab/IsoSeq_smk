# rules/common.smk

import os
import gzip

# Check if "deepvariant_vcf" exists for an individual
def has_vcf(wildcards):
    vcf_path = get_vcf_path(wildcards)
    return bool(vcf_path) and os.path.exists(vcf_path)

# Get modified phased vcf file path
def get_mod_phased_vcf(wildcards):
    phased_vcf = f"mod_vcf/{wildcards.individual}_mod.vcf.gz"
    return phased_vcf if os.path.exists(phased_vcf) else None

# Get raw VCF path for an individual
def get_vcf_path(wildcards):
    individual = wildcards.individual
    try:
        vcf = config["individuals"].get(individual, {}).get("deepvariant_vcf", "")
        if vcf and os.path.exists(vcf):
            return vcf
        else:
            return f"mod_vcf/{individual}_empty.vcf.gz"
    except Exception as e:
        print(f"[WARNING] get_vcf_path failed for {individual}: {e}")
        return f"mod_vcf/{individual}_empty.vcf.gz"


# Get input files for merging per label (safe)
def get_merge_input(wc):
    try:
        return config["individuals"][wc.individual][wc.condition][wc.label]
    except KeyError:
        print(f"[WARNING] Skipping: {wc.individual} / {wc.condition} / {wc.label}")
        return []

# Combine all aligned BAMs for an individual, if conditions/labels exist
def combine_labels(wildcards):
    individual = wildcards.individual
    result = []
    indiv_data = config["individuals"].get(individual, {})
    for condition in ["treated", "untreated"]:
        if condition in indiv_data and isinstance(indiv_data[condition], dict):
            for label in indiv_data[condition]:
                result.append(f"merged/{individual}_{condition}_{label}_labeled.bam")
    return result

# Get all haplotag info files
def get_haplotag_files(wildcards):
    files = []
    for individual, indiv_data in config["individuals"].items():
        for condition in ["treated", "untreated"]:
            if condition in indiv_data and isinstance(indiv_data[condition], dict):
                for label in indiv_data[condition]:
                    files.append(f"whatshap/{individual}_{condition}_{label}.haplotagged.txt")
    return files

# Get whatshap outputs
def whatshap_outs(wc):
    result = []
    for individual, indiv_data in config["individuals"].items():
        for condition in ["treated", "untreated"]:
            if condition in indiv_data and isinstance(indiv_data[condition], dict):
                for label in indiv_data[condition]:
                    result.append(f"whatshap/{individual}_{condition}_{label}.tagged.bam")
    return result
