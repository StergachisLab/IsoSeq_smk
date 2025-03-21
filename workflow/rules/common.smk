import os

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
    vcf = config["individuals"].get(wildcards.individual, {}).get("deepvariant_vcf", "")
    return vcf if os.path.exists(vcf) else ""

# Get input files for merging per label (with protection)
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
        if condition in indiv_data:
            for label in indiv_data[condition]:
                result.append(f"merged/{individual}_{condition}_{label}_labeled.bam")
    return result

# Get all haplotag info files (optional condition fallback)
def get_haplotag_files(wildcards):
    files = []
    for individual, indiv_data in config["individuals"].items():
        for condition in ["treated", "untreated"]:
            for label in indiv_data.get(condition, {}):
                files.append(f"whatshap/{individual}_{condition}_{label}.haplotagged.txt")
    return files

# Get whatshap outputs (tagged BAMs)
def whatshap_outs(wc):
    result = []
    for individual, indiv_data in config["individuals"].items():
        for condition in ["treated", "untreated"]:
            for label in indiv_data.get(condition, {}):
                result.append(f"whatshap/{individual}_{condition}_{label}.tagged.bam")
    return result
