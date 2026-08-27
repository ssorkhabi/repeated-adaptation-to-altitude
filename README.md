# Repeated Adaptation to Altitude in Andean Butterflies

### Abstract
Historical constraints often limit repeatability, with the tension extending to the genetic level. Although some loci repeatedly underlie phenotypic evolution, forming “genetic hotspots” or showing gene reuse, the predictability of evolution at the genetic level varies across taxa. Adaptation to elevational gradients presents a natural laboratory to test this, as distantly related species must contend with similar environmental stresses. Here, I analyse genomic data from thousands of individuals across 11 butterfly species spanning two clades that diverged ~77 million years ago. I quantify genetic repeatability by identifying loci associated with altitude within each species and testing whether the same gene families are repeatedly implicated across independent lineages. My results demonstrate significant statistical evidence for genetic repeatability, identifying 41 orthogroups disproportionately associated with altitude across all 11 species. Additionally, using tissue expression data, I find that orthogroups with stronger evidence for repeatability tend to be expressed in fewer tissue types; however, this relationship is not significant. These results demonstrate that genetic repeatability in altitude adaptation is conserved across diverged lineages. My findings can have implications in understanding Lepidopteran responses to environmental stressors that mimic elevational gradients, such as climate change. **This repository contains all the scripts used in this project.**


### Contents
_**Note.**_ Only the first two directory levels are shown here.

```
.
├── ANGSD
│   ├── angsd-lrt-to-bed.sh
│   ├── angsd-saf-maf-gl.sh
│   ├── angsd_assoc.sh
│   ├── compare_pvals.sh
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
