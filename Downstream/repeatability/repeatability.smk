"""
Altitude GWAS → SNP–gene mapping → PicMin → BLAST
Species: Hera, Hmel, Hnum, Hsar, Hana, Isal, Mlys, Mpol, Mmes, Mmot, Mmen
Traits: continuous, relative (bottom 25% vs top 25% altitude per species), and fixed threshold (200m and 1000m)

Rule 6 (PicMin) and all associated scripts are based on https://github.com/TBooker/PicMin and Booker et al (2023)
"""
import os
import pandas as pd
from pathlib import Path
import sys

configfile: "config/config.yaml"

SPECIES      = config["species"]
OF_SPECIES     = config["full_species"]
ROOT_DIR     = Path(config["root_dir"])
PIPELINE_DIR = Path(config["pipeline_dir"])
GWAS_RES     = Path(config["gwas_results"])
PICMIN_RES   = Path(config["picmin_results"])
SNP_GENE     = Path(config["snp_gene_results"])
OG_FILE      = Path(config["orthogroups"])
PHENO_DIR    = Path(config["phenotype_dir"])
PROT_DIR     = Path(config["proteomes_dir"])
BLAST_RES    = Path(config["blast_results"])
TRAITS       = ["continuous", "categorical", "threshold"]

# Final rule
localrules: all

rule all:
    input:
        expand(
            GWAS_RES / "{sp}_filtered.vcf.gz",
            sp=SPECIES,
        ),
        expand(
            PICMIN_RES / "picmin_{trait}_results.csv",
            trait=TRAITS,
        ),
        expand(
            PICMIN_RES / "picmin_{trait}_plots.done",
            trait=TRAITS,
        ),
        expand(
            BLAST_RES / "{of_sp}_vs_dmel_best_hits.tsv", 
            of_sp=OF_SPECIES
        ),
        BLAST_RES / "summary/species_identity_summary.csv",
        expand(
            BLAST_RES / "reciprocal/{trait}/reciprocal_blast_confidence.csv",
            trait=TRAITS),
        expand(
            GWAS_RES / "{sp}_imputed.vcf.gz",
            sp=SPECIES
        ),


# 1. VCF pre-processing
rule filter_vcf:
    input:
        vcf=ROOT_DIR / "{sp}/09_final_vcf/{sp}_full_concatenated.vcf.gz",
    output:
        vcf=GWAS_RES / "{sp}_filtered.vcf.gz",
        tbi=GWAS_RES / "{sp}_filtered.vcf.gz.tbi",
    conda:
        PIPELINE_DIR / "envs/vcf.yml"
    log:
        "logs/{sp}_filter_vcf.log",
    shell:
        """
        bash scripts/filter_vcf.sh \
            {input.vcf} \
            {output.vcf} \
            {config[min_depth]} \
            {config[max_missing]} \
            {config[min_maf]} \
            > {log} 2>&1
        """


# 1b. Extract scaffold list from filtered VCF 
checkpoint get_scaffolds:
    input:
        vcf = GWAS_RES / "{sp}_filtered.vcf.gz",
    output:
        scaffolds = GWAS_RES / "{sp}_scaffolds.txt",
    log:
        "logs/{sp}_get_scaffolds.log",
    conda:
        PIPELINE_DIR / "envs/vcf.yml",
    shell:
        """
        bcftools index --stats {input.vcf} \
            | awk '$3 >= {config[impute_min_snps]} {{print $1}}' \
            | sed 's/|/__/g' \
            > {output.scaffolds} 2> {log}
        echo "Scaffolds: $(wc -l < {output.scaffolds})" >> {log}
        """


# 1c. Impute a single scaffold
rule impute_scaffold:
    input:
        vcf       = GWAS_RES / "{sp}_filtered.vcf.gz",
        scaffolds = GWAS_RES / "{sp}_scaffolds.txt",
    output:
        vcf = temp(GWAS_RES / "imputed_scaffolds/{sp}/{scaffold}.vcf.gz"),
        tbi = temp(GWAS_RES / "imputed_scaffolds/{sp}/{scaffold}.vcf.gz.tbi"),
    log:
        "logs/{sp}_impute_{scaffold}.log",
    wildcard_constraints:
        scaffold = "[^/]+",
    params:
        out_prefix  = str(GWAS_RES / "imputed_scaffolds/{sp}/{scaffold}"),
        memory_gb   = 108
    conda:
        PIPELINE_DIR / "envs/beagle.yml",
    shell:
        """
        bash scripts/impute_scaffold.sh \
            {wildcards.sp} \
            {wildcards.scaffold} \
            {input.vcf} \
            {params.out_prefix} \
            {threads} \
            {params.memory_gb} \
            > {log} 2>&1
        """


# Helper function to resolve scaffold list after checkpoint 
def get_imputed_scaffolds(wildcards):
    checkpoints.get_scaffolds.get(sp=wildcards.sp)
    scaffold_file = str(GWAS_RES / f"{wildcards.sp}_scaffolds.txt")
    with open(scaffold_file) as f:
        scaffolds = [line.strip() for line in f if line.strip()]
    return [
        str(GWAS_RES / "imputed_scaffolds" / wildcards.sp / f"{s}.vcf.gz")
        for s in scaffolds
    ]


# 1d. Merge all per-scaffold imputed VCFs
rule merge_imputed:
    input:
        vcfs = get_imputed_scaffolds,
    output:
        vcf = GWAS_RES / "{sp}_imputed.vcf.gz",
        tbi = GWAS_RES / "{sp}_imputed.vcf.gz.tbi",
    log:
        "logs/{sp}_merge_imputed.log",
    threads:
        16
    resources:
        mem_mb  = 53920,
        runtime = 120,
    conda:
        PIPELINE_DIR / "envs/vcf.yml",
    shell:
        """
        bcftools concat --naive {input.vcfs} \
            -O z -o {output.vcf} --threads 4 > {log} 2>&1
        bcftools index -t {output.vcf} >> {log} 2>&1
        echo "SNPs: $(bcftools stats {output.vcf} | grep 'number of SNPs' | awk '{{print $NF}}')" >> {log}
        echo "Samples: $(bcftools query -l {output.vcf} | wc -l)" >> {log}
        """


# 2a. Generate categorical phenotype file (per-species IQR)
rule make_categorical_phenotypes:
    input:
        pheno=PHENO_DIR / "{sp}_phenotypes.txt",
    output:
        b2=PHENO_DIR / "{sp}_categorical_phenotypes.txt",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/{sp}_categorical_phenotypes.log",
    shell:
        """
        bash scripts/make_categorical_phenotypes.sh \
            {input.pheno} \
            {output.b2} \
            > {log} 2>&1
        """


# 2b. Generate categorical phenotype file (hardcoded threshold)
rule make_threshold_phenotypes:
    input:
        pheno = PHENO_DIR / "{sp}_phenotypes.txt",
    output:
        pheno = PHENO_DIR / "{sp}_threshold_phenotypes.txt",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/{sp}_make_threshold_phenotypes.log",
    params:
        low_max  = config.get("altitude_low_max",  200),
        high_min = config.get("altitude_high_min", 1000),
    shell:
        """
        bash scripts/make_threshold_phenotypes.sh \
            {input.pheno} \
            {output.pheno} \
            {params.low_max} \
            {params.high_min} \
            > {log} 2>&1
        """

# 3a. PLINK + GEMMA: continuous GWAS
rule gwas_continuous:
    input:
        vcf=GWAS_RES / "{sp}_imputed.vcf.gz",
        pheno=PHENO_DIR / "{sp}_phenotypes.txt",
    output:
        assoc=GWAS_RES / "{sp}_altitude_continuous.assoc.gemma.assoc.txt",
        hits=GWAS_RES  / "{sp}_altitude_continuous_hits.txt",
    log:
        "logs/{sp}_continuous_gwas.log",
    threads: config["threads"]
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    shell:
        """
        bash scripts/gwas_continuous.sh \
            {wildcards.sp} \
            {input.vcf} \
            {input.pheno} \
            {GWAS_RES} \
            {config[n_pcs]} \
            {config[missingness]} \
            {threads} \
            > {log} 2>&1
        """


# 3b. PLINK + GEMMA: categorical GWAS (bottom 25% vs top 25% per species)
rule gwas_categorical:
    input:
        vcf=GWAS_RES / "{sp}_imputed.vcf.gz",
        pheno=PHENO_DIR / "{sp}_categorical_phenotypes.txt",
    output:
        assoc=GWAS_RES / "{sp}_altitude_categorical.assoc.gemma.assoc.txt",
        hits=GWAS_RES  / "{sp}_altitude_categorical_hits.txt",
    log:
        "logs/{sp}_categorical_gwas.log",
    threads: config["threads"]
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    shell:
        """
        bash scripts/gwas_categorical.sh \
            {wildcards.sp} \
            {input.vcf} \
            {input.pheno} \
            {GWAS_RES} \
            {config[n_pcs]} \
            {config[missingness]} \
            {threads} \
            > {log} 2>&1
        """

# 3c. PLINK + GEMMA: harcoded threshold GWAS
rule gwas_threshold:
    input:
        vcf   = GWAS_RES / "{sp}_imputed.vcf.gz",   
        pheno = PHENO_DIR / "{sp}_threshold_phenotypes.txt",
    output:
        assoc = GWAS_RES / "{sp}_altitude_threshold.assoc.gemma.assoc.txt",
        bed   = GWAS_RES / "{sp}_altitude_threshold.bed",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/{sp}_gwas_threshold.log",
    shell:
        """
        bash scripts/gwas_threshold.sh \
            {wildcards.sp} \
            {input.vcf} \
            {input.pheno} \
            {GWAS_RES} \
            {config[n_pcs]} \
            {config[missingness]} \
            {threads} \
            > {log} 2>&1
        """

# 4. Gene BED from GFF (once per species, shared across gwas)
rule make_gene_bed:
    input:
        gff_dir=ROOT_DIR / "{sp}/03_genome",
    output:
        bed=SNP_GENE / "genes_{sp}.bed",
    log:
        "logs/{sp}_make_gene_bed.log",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    shell:
        """
        bash scripts/make_gene_bed.sh \
            {input.gff_dir} \
            {output.bed} \
            {wildcards.sp} \
            "{config[gene_id_strategy_str]}" \
            > {log} 2>&1
        """


# 5. SNP BED + bedtools intersect (per species × trait)
rule snp_gene_map:
    input:
        assoc=GWAS_RES / "{sp}_altitude_{trait}.assoc.gemma.assoc.txt",
        genes=SNP_GENE / "genes_{sp}.bed",
    output:
        snps=SNP_GENE / "snps_{sp}_{trait}.bed",
        map=SNP_GENE  / "snp_gene_map_{sp}_{trait}.txt",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/{sp}_{trait}_snp_gene_map.log"
    shell:
        """
        bash scripts/snp_gene_map.sh \
            {input.assoc} \
            {input.genes} \
            {output.snps} \
            {output.map} \
            > {log} 2>&1
        """

# 6. PicMin (once per trait, across all species)
rule picmin:
    input:
        assocs=expand(
            GWAS_RES / "{sp}_altitude_{{trait}}.assoc.gemma.assoc.txt",
            sp=SPECIES,
        ),
        maps=expand(
            SNP_GENE / "snp_gene_map_{sp}_{{trait}}.txt",
            sp=SPECIES,
        ),
        orthogroups=OG_FILE,
    output:
        results=PICMIN_RES / "picmin_{trait}_results.csv",
        done=PICMIN_RES    / "picmin_{trait}_plots.done",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/picmin_{trait}.log",
    shell:
        """
        bash scripts/picmin.sh \
            {wildcards.trait} \
            {GWAS_RES} \
            {SNP_GENE} \
            {PICMIN_RES} \
            {input.orthogroups} \
            "{config[species_str]}" \
            "{config[orthofinder_names_str]}" \
            {config[picmin_alpha_adapt]} \
            {config[picmin_num_reps]} \
            {config[picmin_null_reps]} \
            {config[picmin_fdr_threshold]} \
            > {log} 2>&1
        """


# 7. Download Drosophila proteome and build BLAST database
rule blast_dmel_db:
    output:
        faa = BLAST_RES / "db/dmel_proteins.faa",
        db  = BLAST_RES / "db/dmel_db.phr",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/blast_dmel_db.log",
    shell:
        """
        bash scripts/blast_dmel_setup.sh {BLAST_RES}/db \
            > {log} 2>&1
        """


# 8. Run blastp for each species
rule blast_dmel_species:
    input:
        faa      = PROT_DIR / "{of_sp}.faa",
        db       = BLAST_RES / "db/dmel_db.phr",
        dmel_faa = BLAST_RES / "db/dmel_proteins.faa",
    output:
        raw  = BLAST_RES / "{of_sp}_vs_dmel_raw.tsv",
        best = BLAST_RES / "{of_sp}_vs_dmel_best_hits.tsv",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/blast_dmel_{of_sp}.log",
    threads:
        8
    resources:
        mem_mb   = 8000,
        runtime  = 60,
    shell:
        """
        bash scripts/blast_dmel_run.sh \
            {wildcards.of_sp} \
            {input.faa} \
            {BLAST_RES}/db \
            {BLAST_RES} \
            {threads} \
            > {log} 2>&1
        """


# 9. Summary comparison across all 11 species
rule blast_dmel_summary:
    input:
        expand(BLAST_RES / "{of_sp}_vs_dmel_best_hits.tsv",
               of_sp=OF_SPECIES),
    output:
        table  = BLAST_RES / "summary/species_identity_summary.csv",
        plot   = BLAST_RES / "summary/pct_identity_by_species.png",
        consist= BLAST_RES / "summary/symbol_consistency.png",
    conda:
        PIPELINE_DIR / "envs/gwas.yml"
    log:
        "logs/blast_dmel_summary.log",
    resources:
        mem_mb  = 4000,
        runtime = 20,
    shell:
        """
        Rscript scripts/blast_dmel_summary.R \
            {BLAST_RES} \
            {BLAST_RES}/summary \
            > {log} 2>&1
        """
