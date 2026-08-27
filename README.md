# Repeated Adaptation to Altitude in Andean Butterflies

### Description
In this project, I analyse genomic data from 1772 individuals **(see `SupplementaryTables`)** across 11 butterfly species spanning two clades that diverged ~77 million years ago. I quantify genetic repeatability by identifying altitude-associated loci within each species and testing whether the same gene families are repeatedly implicated across independent lineages. My results demonstrate significant statistical evidence for genetic repeatability, identifying 41 orthogroups disproportionately associated with altitude across all 11 species **(see `SNP-calling` and `Downstream/repeatability`for scripts)**. Additionally, using tissue expression data, I find that orthogroups with stronger evidence for repeatability tend to be expressed in fewer tissue types; however, this relationship is not significant **(see `Downstream/tissue-specificity` for scripts)**.

For scripts comparing ANGSD and hard-calls, see `ANGSD`. 

### Contents
_**Note.**_ Only the first two directory levels are shown here.

```
.
├── ANGSD
│   ├── angsd-lrt-to-bed.sh
│   ├── angsd-saf-maf-gl.sh
│   ├── angsd_assoc.sh
│   ├── full-comparison.sh
│   ├── gemma-to-bed.sh
│   ├── make-angsd-pheno.sh
│   ├── make-sites-list.sh
│   ├── plot-pval-comparison.R
│   └── plot-pval.sh
├── Downstream
│   ├── repeatability
│   ├── standalone-scripts
│   └── tissue-specificity
├── SNP-calling
│   ├── alignment.smk
│   ├── envs
│   ├── scripts
│   ├── slurm
│   ├── slurm_vcf
│   └── vcf.smk
└── SupplementaryTables
    ├── SupplementaryTable-S2-1.xlsx
    ├── SupplementaryTable-S2-2.xlsx
    └── SupplementaryTable-S2-4.xlsx
```
