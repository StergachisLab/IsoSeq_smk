# rules/common.smk

import os
import gzip
import os, shutil, pathlib, textwrap

#  Node‑local scratch helpers
SCRATCH_PREFIX = os.environ.get("TMPDIR", "/tmp")   # honors scheduler‑set $TMPDIR
PROJECT_ROOT   = workflow.basedir                  # absolute path of repo

def scratch(path):
    """
    Convert a project‑relative path (results/my.bam) to a scratch path
    (/tmp/snakemake/results/my.bam) **without changing directory structure**.
    """
    return os.path.join(SCRATCH_PREFIX, "snakemake", os.path.relpath(path, PROJECT_ROOT))

def stage_to_tmp(files, subdir):
    """
    Symlink every file in *files* into the scratch *subdir* and return the list
    of symlink paths.  Idempotent—re‑running doesn’t recreate links.
    """
    os.makedirs(subdir, exist_ok=True)
    staged = []
    for f in files:
        link = os.path.join(subdir, os.path.basename(f))
        if not os.path.exists(link):
            os.symlink(os.path.abspath(f), link)
        staged.append(link)
    return staged


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
    vcf = config["individuals"].get(individual, {}).get("deepvariant_vcf", "")
    if vcf and os.path.exists(vcf):
        return vcf
    else:
        # Make sure rule `create_dummy_vcf` has this as its output
        return f"mod_vcf/{individual}_empty.vcf.gz"


# Get input files for merging per label (safe)
def get_merge_input(wildcards):
    raw_files = config["individuals"][wildcards.individual][wildcards.condition][wildcards.label]
    return stage_to_tmp(
        raw_files,
        os.path.join(SCRATCH_PREFIX,
                     "snakemake",
                     "raw",
                     wildcards.individual,
                     wildcards.condition,
                     wildcards.label)
    )


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
                    files.append(
                        f"whatshap/{individual}_{condition}_{label}.haplotagged.txt"
                    )
    return files


# Get whatshap outputs
def whatshap_outs(wc):
    result = []
    for individual, indiv_data in config["individuals"].items():
        for condition in ["treated", "untreated"]:
            if condition in indiv_data and isinstance(indiv_data[condition], dict):
                for label in indiv_data[condition]:
                    result.append(
                        f"whatshap/{individual}_{condition}_{label}.tagged.bam"
                    )
    return result
