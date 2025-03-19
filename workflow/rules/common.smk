# check if "deepvariant_vcf" exists for an individual
def has_vcf(wildcards):
    vcf_path = get_vcf_path(wildcards)
    return vcf_path is not None and os.path.exists(vcf_path)


# get modified phased vcf file
def get_mod_phased_vcf(wildcards):
    phased_vcf = f"mod_vcf/{wildcards.individual}_mod.vcf.gz"
    if os.path.exists(phased_vcf):
        return phased_vcf
    return None


# get vcf path for an individual
def get_vcf_path(wildcards):
    vcf = config["individuals"][wildcards["individual"]].get("deepvariant_vcf")
    if vcf and os.path.exists(vcf):
        return vcf
    return []


# get the merged input for an individual, condition, and label
def get_merge_input(wc):
    if wc.label in config["individuals"][wc.individual][wc.condition]:
        return config["individuals"][wc.individual][wc.condition][wc.label]
    else:
        raise KeyError(
            f"Label {wc.label} not found for individual {wc.individual} "
            f"and condition {wc.condition}."
        )


# combine labels for an individual
def combine_labels(wildcards):
    individual = wildcards.individual
    rtn = []
    for condition in ["treated", "untreated"]:
        # for condition in config["individuals"][individual]:
        for label in config["individuals"][individual][condition]:
            file = f"merged/{individual}_{condition}_{label}_labeled.bam"
            rtn.append(file)
    return rtn


# get haplotag files
def get_haplotag_files(wildcards):
    files = []
    for individual in config["individuals"]:
        for condition in ["treated", "untreated"]:
            # for condition in config["individuals"][individual]:
            for label in config["individuals"][individual][condition]:
                files.append(
                    f"whatshap/{individual}_{condition}_{label}.haplotagged.txt"
                )
    return files


# get whatshap output files
def whatshap_outs(wc):
    template_f_path = "whatshap/{individual}_{condition}_{label}.tagged.bam"
    rtn = []
    for individual in config["individuals"]:
        # for condition in config["individuals"][individual]:
        for condition in ["treated", "untreated"]:
            for label in config["individuals"][individual][condition]:
                rtn.append(
                    template_f_path.format(
                        individual=individual, condition=condition, label=label
                    )
                )
    return rtn
