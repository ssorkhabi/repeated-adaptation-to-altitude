"""
Tissue-specificity Pipeline — TPM → tau → orthogroup-level ep-values → Z-scores
Species: Mechanitis messenoides (Mmes) RNA-seq, applied to all orthogroups

Steps:
  1. Salmon quantification (array job per sample)
  2. tximport + tau calculation + orthogroup-level ep-values + Z-scores (R)
  3. Comparison with PicMin results (R)
"""

import os
import pandas as pd
from pathlib import Path
import sys

configfile: "config/config.yaml"

ROOT_DIR     = Path(config["root_dir"])
PIPELINE_DIR = Path(config["pipeline_dir"])
QUANTS_DIR   = Path(config["quants_dir"])
TAU_DIR      = Path(config["tau_dir"])
PLOTS_DIR     = Path(config["plots_dir"])
SAMPLE_LIST  = Path(config["sample_list"])
INDEX        = Path(config["salmon_index"])
GFF          = Path(config["gff"])
METADATA     = config["sample_metadata"]
ORTHOGROUPS  = config["orthogroups"]
PICMIN_CONT  = config["picmin_continuous"]
PICMIN_BIN   = config["picmin_categorical"]


# Read sample list 
with open(SAMPLE_LIST) as f:
    SAMPLES = [s.strip() for s in f if s.strip()]


# Final
rule all:
    input:
        # Salmon quants for all samples
        expand(QUANTS_DIR / "{sample}/quant.sf", sample=SAMPLES),
        # Gene-level tau
        TAU_DIR / "geneLevel_tau.csv",
        TAU_DIR / "mean_tpm.csv",
        # Orthogroup-level ep-values and Z-scores
        TAU_DIR / "orthogroup_tau_Zscore.csv",
        # Diagnostic plots
        PLOTS_DIR / "tau_distribution.png",
        PLOTS_DIR / "TPM_by_tissue.png",
        PLOTS_DIR / "orthogroup_epvalue_distribution.png",
        # PicMin comparison
        TAU_DIR / "picmin_pleiotropy_comparison.csv",
        PLOTS_DIR / "picmin_pleiotropy_boxplot.png",
        # Figues
        PLOTS_DIR / "tissue_expression_significant_all11.png",
        PLOTS_DIR / "coexpression_network_centrality.png",
        PLOTS_DIR / "pleiotropy_duplication_barplot.png",


# 1. Salmon quantification (one job per sample)
rule salmon_quant:
    input:
        r1=PIPELINE_DIR / "rawData/{sample}_1.fastq.gz",
        r2=PIPELINE_DIR / "rawData/{sample}_2.fastq.gz",
        index=INDEX,
    output:
        quant=QUANTS_DIR / "{sample}/quant.sf",
    log:
        "logs/salmon_{sample}.log",
    threads: config["salmon_threads"]
    conda:
        PIPELINE_DIR / "envs/pleiotropy.yml"
    shell:
        """
        bash scripts/salmon_quant.sh \
            {wildcards.sample} \
            {input.r1} \
            {input.r2} \
            {input.index} \
            {QUANTS_DIR} \
            {threads} \
            > {log} 2>&1
        """


# 2. TPM → tau → orthogroup ep-values + Z-scores
rule compute_tau:
    input:
        quants=expand(QUANTS_DIR / "{sample}/quant.sf", sample=SAMPLES),
        gff=GFF,
        metadata=METADATA,
        orthogroups=ORTHOGROUPS,
    output:
        mean_tpm=TAU_DIR / "mean_tpm.csv",
        gene_tau=TAU_DIR / "geneLevel_tau.csv",
        og_tau=TAU_DIR / "orthogroup_tau_Zscore.csv",
        plot_tau=PLOTS_DIR / "tau_distribution.png",
        plot_tpm=PLOTS_DIR / "TPM_by_tissue.png",
        plot_ep=PLOTS_DIR / "orthogroup_epvalue_distribution.png",
    log:
        "logs/compute_tau.log",
    conda:
        PIPELINE_DIR / "envs/pleiotropy.yml"
    shell:
        """
        bash scripts/compute_tau.sh \
            {QUANTS_DIR} \
            {input.gff} \
            {input.metadata} \
            {input.orthogroups} \
            {TAU_DIR} \
            {PLOTS_DIR} \
            {config[of_col_mmes]} \
            > {log} 2>&1
        """


# 3. Compare orthogroup tau Z-scores with PicMin results
rule picmin_pleiotropy_comparison:
    input:
        og_tau=TAU_DIR / "orthogroup_tau_Zscore.csv",
        picmin_cont=PICMIN_CONT,
        picmin_bin=PICMIN_BIN,
    output:
        results=TAU_DIR / "picmin_pleiotropy_comparison.csv",
        plot=PLOTS_DIR / "picmin_pleiotropy_boxplot.png",
    log:
        "logs/picmin_pleiotropy_comparison.log",
    conda:
        PIPELINE_DIR / "envs/pleiotropy.yml"
    shell:
        """
        bash scripts/picmin_pleiotropy_comparison.sh \
            {input.og_tau} \
            {input.picmin_cont} \
            {input.picmin_bin} \
            {TAU_DIR} \
            {PLOTS_DIR} \
            {config[picmin_fdr_threshold]} \
            > {log} 2>&1
        """


# 4. Tissue expression panel for orthogroups significant in all 11 species
rule tissue_expression_significant:
    input:
        og_tau=TAU_DIR / "orthogroup_tau_Zscore.csv",
        mean_tpm=TAU_DIR / "mean_tpm.csv",
        gene_tau=TAU_DIR / "geneLevel_tau.csv",
        picmin_cont=PICMIN_CONT,
        picmin_bin=PICMIN_BIN,
    output:
        plot=PLOTS_DIR / "tissue_expression_significant_all11.png",
        summary=TAU_DIR / "tissue_expression_significant_all11_summary.csv",
    log:
        "logs/tissue_expression_significant.log",
    conda:
        PIPELINE_DIR / "envs/pleiotropy.yml"
    params:
        fdr=config["picmin_fdr_threshold"],
        max_ogs=config.get("max_significant_ogs_plotted", 20),
        picmin_dir=lambda wc: PICMIN_CONT.rsplit("/", 1)[0],
    shell:
        """
        cd {workflow.basedir} && \
        PICMIN_DIR={params.picmin_dir} \
        ORTHOGROUPS_PATH={ORTHOGROUPS} \
        OF_COL_MMES={config[of_col_mmes]} \
        TAU_DIR_OVERRIDE={TAU_DIR} \
        Rscript scripts/plot_tissue_expression_significant.R \
            {params.fdr} {params.max_ogs} \
            > {log} 2>&1
        """

