#!/bin/bash

# compare_wald_angsd.sh
#
# Compare Wald/GEMMA and ANGSD GWAS results:
#   1. Exact significant SNP overlap
#   2. -log10(p) regression
#   3. Significant SNP BED files
#   4. SNP -> gene mapping
#   5. Candidate-gene overlap

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate gwas

set -euo pipefail

THREADS=16
TMP="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/tmp"

mkdir -p comparison/results
mkdir -p comparison/bed
mkdir -p comparison/genes

# INPUT FILES
MMOT_FILE="comparison/Mmot_wald_vs_angsd_shared.tsv"
HERA_FILE="comparison/Hera_wald_vs_angsd_shared.tsv"

# BONFERRONI THRESHOLDS
# Wald / GEMMA
MMOT_WALD_THR="1.1056641182583e-09"
HERA_WALD_THR="8.16394414640878e-10"

# ANGSD
# Mmot: 0.05 / 59,030,459
# Hera: 0.05 / 54,552,912
MMOT_ANGSD_THR="8.470203e-10"
HERA_ANGSD_THR="9.165433e-10"

# GENE BED FILES
MMOT_GENE_BED="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/gwas_pipeline/results/snp_gene_map/genes_Mmot.bed"
HERA_GENE_BED="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/gwas_pipeline/results/snp_gene_map/genes_Hera.bed"


set -euo pipefail


mkdir -p comparison/results
mkdir -p comparison/bed
mkdir -p comparison/genes
mkdir -p comparison/support



# INPUT FILES


MMOT_FILE="comparison/Mmot_wald_vs_angsd_shared.tsv"
HERA_FILE="comparison/Hera_wald_vs_angsd_shared.tsv"



# BONFERRONI THRESHOLDS


# Wald / GEMMA
MMOT_WALD_THR="1.1056641182583e-09"
HERA_WALD_THR="8.16394414640878e-10"

# ANGSD
# Mmot: 0.05 / 59,030,459
# Hera: 0.05 / 54,552,912
MMOT_ANGSD_THR="8.470203e-10"
HERA_ANGSD_THR="9.165433e-10"



# GENE BED FILES
# EDIT THESE
#
# Expected format:
# chrom   start   end   gene_id


MMOT_GENE_BED="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/gwas_pipeline/results/snp_gene_map/genes_Mmot.bed"
HERA_GENE_BED="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/gwas_pipeline/results/snp_gene_map/genes_Hera.bed"


# TEMP DIRECTORY
TMP="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/tmp"

mkdir -p "$TMP"



# PROCESS EACH SPECIES


process_species () {

    NAME=$1
    FILE=$2
    WALD_THR=$3
    ANGSD_THR=$4
    GENE_BED=$5

    echo "Processing $NAME"


    
    # 1. Significant Wald SNPs -> BED
    

    awk -v p="$WALD_THR" '
        $3 < p {
            print $1, $2-1, $2, $1 ":" $2, $3
        }
    ' OFS='\t' "$FILE" \
    | sort -T "$TMP" -k1,1 -k2,2n \
    > "comparison/bed/${NAME}_wald_sig.bed"


    
    # 2. Significant ANGSD SNPs -> BED
    

    awk -v p="$ANGSD_THR" '
        $4 < p {
            print $1, $2-1, $2, $1 ":" $2, $4
        }
    ' OFS='\t' "$FILE" \
    | sort -T "$TMP" -k1,1 -k2,2n \
    > "comparison/bed/${NAME}_angsd_sig.bed"


    NW_SNP=$(wc -l < "comparison/bed/${NAME}_wald_sig.bed")
    NA_SNP=$(wc -l < "comparison/bed/${NAME}_angsd_sig.bed")

    echo "Wald significant SNPs:  $NW_SNP"
    echo "ANGSD significant SNPs: $NA_SNP"


    
    # 3. Map significant Wald SNPs to genes
    #
    # A has 5 columns.
    # B gene BED has 4 columns.
    # Gene ID therefore = column 9.
    

    bedtools intersect \
        -a "comparison/bed/${NAME}_wald_sig.bed" \
        -b "$GENE_BED" \
        -wa -wb \
    | awk '{print $4, $9}' OFS='\t' \
    > "comparison/genes/${NAME}_wald_snp_gene_map.tsv"


    
    # 4. Map significant ANGSD SNPs to genes
    

    bedtools intersect \
        -a "comparison/bed/${NAME}_angsd_sig.bed" \
        -b "$GENE_BED" \
        -wa -wb \
    | awk '{print $4, $9}' OFS='\t' \
    > "comparison/genes/${NAME}_angsd_snp_gene_map.tsv"


    
    # 5. Unique candidate genes
    

    cut -f2 "comparison/genes/${NAME}_wald_snp_gene_map.tsv" \
        | sort -u \
        > "comparison/genes/${NAME}_wald_genes.txt"

    cut -f2 "comparison/genes/${NAME}_angsd_snp_gene_map.tsv" \
        | sort -u \
        > "comparison/genes/${NAME}_angsd_genes.txt"


    
    # 6. Shared genome-wide significant genes
    

    comm -12 \
        "comparison/genes/${NAME}_wald_genes.txt" \
        "comparison/genes/${NAME}_angsd_genes.txt" \
        > "comparison/genes/${NAME}_shared_genes.txt"


    NW_GENE=$(wc -l < "comparison/genes/${NAME}_wald_genes.txt")
    NA_GENE=$(wc -l < "comparison/genes/${NAME}_angsd_genes.txt")
    NB_GENE=$(wc -l < "comparison/genes/${NAME}_shared_genes.txt")

    echo "Wald candidate genes:   $NW_GENE"
    echo "ANGSD candidate genes:  $NA_GENE"
    echo "Shared candidate genes: $NB_GENE"


    
    # 7. Create BED of ALL shared SNPs carrying ANGSD p-value
    #
    # This is important:
    # we now look at ALL ANGSD p-values inside Wald genes,
    # not only genome-wide significant ANGSD SNPs.
    

    awk '
        {
            print $1, $2-1, $2, $4
        }
    ' OFS='\t' "$FILE" \
    | sort -T "$TMP" -k1,1 -k2,2n \
    > "comparison/bed/${NAME}_all_angsd_p.bed"


    
    # 8. Extract gene coordinates for Wald candidate genes
    

    awk '
        NR==FNR {
            keep[$1]=1
            next
        }

        $4 in keep {
            print
        }
    ' \
        "comparison/genes/${NAME}_wald_genes.txt" \
        "$GENE_BED" \
    > "comparison/genes/${NAME}_wald_candidate_genes.bed"


    
    # 9. Intersect ALL ANGSD SNPs with Wald candidate genes
    #
    # Columns:
    # 1 chrom
    # 2 start
    # 3 end
    # 4 ANGSD p
    # 5 gene chrom
    # 6 gene start
    # 7 gene end
    # 8 gene ID
    

    bedtools intersect \
        -a "comparison/bed/${NAME}_all_angsd_p.bed" \
        -b "comparison/genes/${NAME}_wald_candidate_genes.bed" \
        -wa -wb \
    > "comparison/support/${NAME}_ANGSD_in_Wald_genes.tsv"


    
    # 10. Gene overlap summary
    

    {
        echo -e "Species\t$NAME"
        echo -e "Wald_threshold\t$WALD_THR"
        echo -e "ANGSD_threshold\t$ANGSD_THR"

        echo -e "Wald_significant_SNPs\t$NW_SNP"
        echo -e "ANGSD_significant_SNPs\t$NA_SNP"

        echo -e "Wald_candidate_genes\t$NW_GENE"
        echo -e "ANGSD_candidate_genes\t$NA_GENE"
        echo -e "Shared_candidate_genes\t$NB_GENE"

        awk -v n="$NB_GENE" -v d="$NW_GENE" 'BEGIN {
            if (d > 0)
                print "Proportion_Wald_genes_shared\t" n/d
            else
                print "Proportion_Wald_genes_shared\tNA"
        }'

        awk -v n="$NB_GENE" -v d="$NA_GENE" 'BEGIN {
            if (d > 0)
                print "Proportion_ANGSD_genes_shared\t" n/d
            else
                print "Proportion_ANGSD_genes_shared\tNA"
        }'

    } > "comparison/results/${NAME}_gene_overlap.txt"
}



# RUN BOTH SPECIES


process_species \
    "Mmot" \
    "$MMOT_FILE" \
    "$MMOT_WALD_THR" \
    "$MMOT_ANGSD_THR" \
    "$MMOT_GENE_BED"

process_species \
    "Hera" \
    "$HERA_FILE" \
    "$HERA_WALD_THR" \
    "$HERA_ANGSD_THR" \
    "$HERA_GENE_BED"



# R ANALYSIS


Rscript - <<EOF

library(ggplot2)

datasets <- list(

    Mmot = list(
        file = "$MMOT_FILE",
        wald_thr = as.numeric("$MMOT_WALD_THR"),
        angsd_thr = as.numeric("$MMOT_ANGSD_THR")
    ),

    Hera = list(
        file = "$HERA_FILE",
        wald_thr = as.numeric("$HERA_WALD_THR"),
        angsd_thr = as.numeric("$HERA_ANGSD_THR")
    )
)


for (name in names(datasets)) {

    z <- datasets[[name]]

    d <- read.table(z\$file)


    
    # 1. Genome-wide p-value comparison
    

    d\$wald_logp  <- -log10(d\$V3)
    d\$angsd_logp <- -log10(d\$V4)

    d_plot <- d[
        is.finite(d\$wald_logp) &
        is.finite(d\$angsd_logp),
    ]


    p <- ggplot(
        d_plot,
        aes(wald_logp, angsd_logp)
    ) +
        geom_point(
            alpha = 0.3,
            size = 0.5
        ) +
        geom_smooth(
            method = "lm",
            se = FALSE,
            color = "red"
        ) +
        theme_classic() +
        labs(
            x = "-log10(Wald p)",
            y = "-log10(ANGSD p)",
            title = name
        )


    ggsave(
        paste0(
            "comparison/results/",
            name,
            "_pvalue_comparison.png"
        ),
        p,
        width = 5,
        height = 5,
        dpi = 300
    )


    
    # 2. Regression
    

    fit <- lm(
        angsd_logp ~ wald_logp,
        data = d_plot
    )

    capture.output(
        summary(fit),
        file = paste0(
            "comparison/results/",
            name,
            "_regression.txt"
        )
    )


    
    # 3. Exact genome-wide significant SNP overlap
    

    wald_sig  <- d\$V3 < z\$wald_thr
    angsd_sig <- d\$V4 < z\$angsd_thr

    n_wald <- sum(
        wald_sig,
        na.rm = TRUE
    )

    n_angsd <- sum(
        angsd_sig,
        na.rm = TRUE
    )

    n_both <- sum(
        wald_sig & angsd_sig,
        na.rm = TRUE
    )


    prop_wald <- if (n_wald > 0) {
        n_both / n_wald
    } else {
        NA
    }

    prop_angsd <- if (n_angsd > 0) {
        n_both / n_angsd
    } else {
        NA
    }


    overlap_table <- table(
        Wald = wald_sig,
        ANGSD = angsd_sig
    )


    capture.output({

        cat("Species:", name, "\n\n")

        cat(
            "Wald Bonferroni threshold:",
            z\$wald_thr,
            "\n"
        )

        cat(
            "ANGSD Bonferroni threshold:",
            z\$angsd_thr,
            "\n\n"
        )

        cat(
            "Significant Wald SNPs:",
            n_wald,
            "\n"
        )

        cat(
            "Significant ANGSD SNPs:",
            n_angsd,
            "\n"
        )

        cat(
            "Significant in both:",
            n_both,
            "\n\n"
        )

        cat(
            "Proportion Wald -> ANGSD:",
            prop_wald,
            "\n"
        )

        cat(
            "Proportion ANGSD -> Wald:",
            prop_angsd,
            "\n\n"
        )

        cat("Contingency table:\n")

        print(overlap_table)

    },
    file = paste0(
        "comparison/results/",
        name,
        "_SNP_overlap.txt"
    ))


    
    # 4. ANGSD support within Wald candidate genes
    

    support_file <- paste0(
        "comparison/support/",
        name,
        "_ANGSD_in_Wald_genes.tsv"
    )


    if (file.exists(support_file) &&
        file.info(support_file)\$size > 0) {

        x <- read.table(
            support_file,
            header = FALSE
        )


        colnames(x) <- c(
            "chrom",
            "start",
            "end",
            "angsd_p",
            "gene_chr",
            "gene_start",
            "gene_end",
            "gene"
        )


        
        # Minimum ANGSD p-value for each Wald candidate gene
        

        gene_support <- aggregate(
            angsd_p ~ gene,
            data = x,
            FUN = min
        )


        
        # Number of ANGSD SNPs tested in each gene
        

        n_snps <- aggregate(
            angsd_p ~ gene,
            data = x,
            FUN = length
        )

        colnames(n_snps)[2] <- "n_angsd_snps"


        gene_support <- merge(
            gene_support,
            n_snps,
            by = "gene"
        )


        colnames(gene_support)[
            colnames(gene_support) == "angsd_p"
        ] <- "min_angsd_p"


        
        # Nominal support
        

        gene_support\$ANGSD_p_lt_0.05 <-
            gene_support\$min_angsd_p < 0.05

        gene_support\$ANGSD_p_lt_0.01 <-
            gene_support\$min_angsd_p < 0.01

        gene_support\$ANGSD_p_lt_0.001 <-
            gene_support\$min_angsd_p < 0.001


        
        # Genome-wide ANGSD significance
        

        gene_support\$ANGSD_genomewide_sig <-
            gene_support\$min_angsd_p < z\$angsd_thr


        
        # Correct across Wald candidate genes
        #
        # These are exploratory locus-level summaries.
        

        gene_support\$BH_gene_p <-
            p.adjust(
                gene_support\$min_angsd_p,
                method = "BH"
            )

        gene_support\$Bonferroni_gene_p <-
            p.adjust(
                gene_support\$min_angsd_p,
                method = "bonferroni"
            )


        gene_support <- gene_support[
            order(gene_support\$min_angsd_p),
        ]


        
        # Save per-gene table
        

        write.table(
            gene_support,
            paste0(
                "comparison/results/",
                name,
                "_Wald_genes_ANGSD_support.tsv"
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )


        
        # Summary
        

        n_genes <- nrow(gene_support)

        n_005 <- sum(
            gene_support\$ANGSD_p_lt_0.05
        )

        n_001 <- sum(
            gene_support\$ANGSD_p_lt_0.01
        )

        n_0001 <- sum(
            gene_support\$ANGSD_p_lt_0.001
        )

        n_genome <- sum(
            gene_support\$ANGSD_genomewide_sig
        )

        n_bh <- sum(
            gene_support\$BH_gene_p < 0.05
        )

        n_bonf_gene <- sum(
            gene_support\$Bonferroni_gene_p < 0.05
        )


        capture.output({

            cat(
                "Species:",
                name,
                "\n\n"
            )

            cat(
                "Wald candidate genes with ANGSD SNPs:",
                n_genes,
                "\n\n"
            )

            cat(
                "Genes with min ANGSD p < 0.05:",
                n_005,
                "(",
                n_005 / n_genes,
                ")\n"
            )

            cat(
                "Genes with min ANGSD p < 0.01:",
                n_001,
                "(",
                n_001 / n_genes,
                ")\n"
            )

            cat(
                "Genes with min ANGSD p < 0.001:",
                n_0001,
                "(",
                n_0001 / n_genes,
                ")\n\n"
            )

            cat(
                "Genes containing a genome-wide significant ANGSD SNP:",
                n_genome,
                "(",
                n_genome / n_genes,
                ")\n\n"
            )

            cat(
                "Genes significant after BH correction across Wald candidate genes:",
                n_bh,
                "(",
                n_bh / n_genes,
                ")\n"
            )

            cat(
                "Genes significant after Bonferroni correction across Wald candidate genes:",
                n_bonf_gene,
                "(",
                n_bonf_gene / n_genes,
                ")\n"
            )

        },
        file = paste0(
            "comparison/results/",
            name,
            "_Wald_genes_ANGSD_support_summary.txt"
        ))
    }
}

EOF

echo "DONEEE!"