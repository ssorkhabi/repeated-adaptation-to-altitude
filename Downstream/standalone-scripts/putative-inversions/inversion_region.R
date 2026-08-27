#!/usr/bin/env Rscript
# inversion_region.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(ggplot2)
  library(ggtext)
  library(snpStats)
  library(LDheatmap)
  library(cowplot)
  library(magick)
})

# 0. Args
valid_sp <- c("Hera", "Hsar", "Hmel", "Hnum", "Mlys", "Mpol", "Mmes", "Mmot", "Mmen", "Isal", "Hana")

## Species colours
this_script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(this_script_dir) == 0 || this_script_dir == "") this_script_dir <- "."
source(file.path(this_script_dir, "species_config.R"))

option_list <- list(
  make_option("--sp",           type = "character", help = "Species code (Hana, Hera, Hnum, Mmes)"),
  make_option("--scaffold",     type = "character", default = NA_character_, help = "Target raw scaffold name for a zoomed LD/Manhattan region. Mutually exclusive with --chr-label."),
  make_option("--chr-label",    type = "character", default = NA_character_, dest = "chr_label",
              help = "Plotted chromosome label (e.g. '2'), using --lookup. With --start/--end omitted: whole-chromosome PREVIEW plot only, no LD heatmap. With --start/--end given (in CONCATENATED chromosome-group coordinates, i.e. raw position + that scaffold's offset from the lookup table): full LD heatmap/Manhattan/composite pipeline spanning potentially multiple scaffolds. Mutually exclusive with --scaffold."),
  make_option("--lookup",       type = "character", default = NA_character_, help = "Path to the _chr_lookup.csv from manhattan.R. Required with --chr-label."),
  make_option("--start",        type = "integer",   default = NA_integer_, help = "Region start (bp). With --scaffold: native scaffold coordinate. With --chr-label: CONCATENATED chromosome-group coordinate (raw position + offset)."),
  make_option("--end",          type = "integer",   default = NA_integer_, help = "Region end (bp). Same coordinate system rules as --start."),
  make_option("--trait",        type = "character", default = "altitude", help = "Trait label for plot titles [default %default]"),
  make_option("--assoc",        type = "character", default = NA_character_, help = "Path to GEMMA assoc file. Default: results/gwas/<sp>_altitude_categorical.assoc.gemma.assoc.txt"),
  make_option("--plink-prefix", type = "character", default = NA_character_, dest = "plink_prefix", help = "PLINK bed/bim/fam prefix (no extension). Default: results/gwas/<sp>_altitude_categorical"),
  make_option("--outdir",       type = "character", default = "results/putative-inversions/inversion-regions", help = "Output directory [default %default]"),
  make_option("--pval-col",     type = "character", default = "p_wald", dest = "pval_col", help = "P-value column in assoc file [default %default]"),
  make_option("--max-ld-snps",  type = "integer",   default = 250, dest = "max_ld_snps",
              help = "After QC, if the window contains more SNPs than this, evenly thin by physical position before computing LD [default %default]"),
  make_option("--min-maf",      type = "double",    default = 0.05, dest = "min_maf",
              help = "Minimum minor allele frequency for SNPs used in the LD heatmap [default %default]"),
  make_option("--min-call-rate", type = "double",   default = 0.90, dest = "min_call_rate",
              help = "Minimum SNP call rate for SNPs used in the LD heatmap [default %default]"),
  make_option("--ld-distance-bin-kb", type = "double", default = 50, dest = "ld_distance_bin_kb",
              help = "Bin width (kb) for the LD-decay plot summarising R^2 by physical distance [default %default]"),
  make_option("--mark-start",   type = "integer",   default = NA_integer_, dest = "mark_start",
              help = "Optional: bp position of a putative inversion breakpoint to flag on the LD heatmap (same coordinate system as --start/--end). Nearest actual marker is labelled."),
  make_option("--mark-end",     type = "integer",   default = NA_integer_, dest = "mark_end",
              help = "Optional: bp position of the other putative inversion breakpoint to flag on the LD heatmap.")
)

opt <- parse_args(OptionParser(option_list = option_list))

has_window <- !is.na(opt$start) && !is.na(opt$end)
chr_label_mode <- !is.na(opt$chr_label)
preview_mode <- chr_label_mode && !has_window        # whole-chromosome scan, no LD heatmap
multiregion_mode <- chr_label_mode && has_window     # chr-label + window -> full pipeline, possibly multi-scaffold
single_mode <- !chr_label_mode                       # --scaffold based, existing behaviour

if (is.null(opt$sp)) stop("--sp is required. See header of this script for usage.")
if (!opt$sp %in% valid_sp) stop("--sp must be one of: ", paste(valid_sp, collapse = ", "))

if (chr_label_mode) {
  if (is.na(opt$lookup)) stop("--lookup is required when using --chr-label.")
  if (!file.exists(opt$lookup)) stop("Lookup file not found: ", opt$lookup)
  if (!is.na(opt$scaffold)) stop("--scaffold and --chr-label are mutually exclusive.")
  if (has_window && opt$end <= opt$start) stop("--end must be greater than --start")
} else {
  if (is.na(opt$scaffold) || is.na(opt$start) || is.na(opt$end)) {
    stop("--scaffold, --start and --end are all required (unless using --chr-label).")
  }
  if (opt$end <= opt$start) stop("--end must be greater than --start")
}

assoc_file  <- if (is.na(opt$assoc))        file.path("results/gwas", paste0(opt$sp, "_altitude_categorical.assoc.gemma.assoc.txt")) else opt$assoc
plink_prefix<- if (is.na(opt$plink_prefix)) file.path("results/gwas", paste0(opt$sp, "_altitude_categorical"))          else opt$plink_prefix

bed <- paste0(plink_prefix, ".bed")
bim <- paste0(plink_prefix, ".bim")
fam <- paste0(plink_prefix, ".fam")

required_files <- if (preview_mode) c(assoc_file) else c(assoc_file, bed, bim, fam)
for (f in required_files) {
  if (!file.exists(f)) stop("Required input not found: ", f,
                             "\n(Check you're running from the gwas_pipeline root, or pass --assoc / --plink-prefix explicitly.)")
}
 
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
tag <- case_when(
  preview_mode     ~ sprintf("%s_chr%s_preview", opt$sp, opt$chr_label),
  multiregion_mode ~ sprintf("%s_chr%s_%d-%d", opt$sp, opt$chr_label, opt$start, opt$end),
  TRUE             ~ sprintf("%s_%s_%d-%d", opt$sp, opt$scaffold, opt$start, opt$end)
)

message("Region: ", tag)

sp_col <- species_colors[[opt$sp]]

# 1. Regional GWAS data
gwas <- read.delim(assoc_file, sep = "\t", header = TRUE)
gwas <- gwas[!grepl("^##", gwas$chr), ]
gwas <- gwas[!is.na(gwas[[opt$pval_col]]), ]
gwas$P <- gwas[[opt$pval_col]]
 
n_total_tests <- nrow(gwas)
bonferroni_threshold <- 0.05 / n_total_tests
message("Total SNPs tested genome-wide: ", n_total_tests, " | Bonferroni threshold: ", signif(bonferroni_threshold, 3))
 
if (chr_label_mode) {
  ## Whole chromosome-GROUP: pull every raw scaffold belonging to this
  ## CHR_label from the lookup table, and use its `offset` (within-chromosome,
  ## already computed by manhattan.R) to concatenate them in order.
  lookup <- read.csv(opt$lookup, stringsAsFactors = FALSE)
  group_scaffolds <- lookup %>% filter(CHR_label == opt$chr_label)
  if (nrow(group_scaffolds) == 0) stop("No scaffolds found for CHR_label '", opt$chr_label, "' in ", opt$lookup)
  message("Chromosome '", opt$chr_label, "' is made of ", nrow(group_scaffolds), " scaffold(s): ",
          paste(group_scaffolds$scaffold, collapse = ", "))

  chr_gwas <- gwas %>%
    inner_join(group_scaffolds %>% select(scaffold, offset), by = c("chr" = "scaffold")) %>%
    mutate(ps_adj = ps + offset) %>%
    filter(P > 0) %>%
    arrange(ps_adj)
  if (nrow(chr_gwas) == 0) stop("No SNPs found for any scaffold in chromosome '", opt$chr_label, "'.")
  message("SNPs across whole chromosome '", opt$chr_label, "': ", nrow(chr_gwas))
 
  if (multiregion_mode) {
    region_gwas <- chr_gwas %>% filter(ps_adj >= opt$start, ps_adj <= opt$end)
    if (nrow(region_gwas) == 0) {
      stop("No SNPs found for chromosome '", opt$chr_label, "' in concatenated range [", opt$start, ", ", opt$end, "].")
    }
    region_scaffolds <- sort(unique(region_gwas$chr))
    message("SNPs in LD/highlight window: ", nrow(region_gwas), " spanning scaffold(s): ",
            paste(region_scaffolds, collapse = ", "))
  }
} else {
  ## Whole-chromosome subset -- used for the Manhattan panel
  chr_gwas <- gwas %>%
    filter(chr == opt$scaffold, P > 0) %>%
    arrange(ps) %>%
    mutate(ps_adj = ps)  # no offset needed, single raw scaffold
  if (nrow(chr_gwas) == 0) stop("No SNPs found for scaffold '", opt$scaffold, "' at all. Check the scaffold name matches the 'chr' column exactly.")
  message("SNPs on whole scaffold: ", nrow(chr_gwas))
 
  ## LD-window subset
  region_gwas <- chr_gwas %>% filter(ps >= opt$start, ps <= opt$end)
  if (nrow(region_gwas) == 0) {
    stop("No SNPs found for scaffold '", opt$scaffold, "' in [", opt$start, ", ", opt$end,
         "]. Check the scaffold name matches the 'chr' column of the assoc file exactly.")
  }
  message("SNPs in LD/highlight window: ", nrow(region_gwas))
}

# 2. Whole-chromosome Manhattan panel
manhattan_png <- file.path(opt$outdir, paste0("manhattan_", tag, ".png"))
ymax <- ceiling(max(-log10(chr_gwas$P), na.rm = TRUE)) + 1

p_manhattan <- ggplot(chr_gwas, aes(x = ps_adj, y = -log10(P)))

if (!preview_mode) {
  p_manhattan <- p_manhattan +
    geom_rect(data = data.frame(xmin = opt$start, xmax = opt$end, ymin = 0, ymax = ymax),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = sp_col, alpha = 0.15)
}

p_manhattan <- p_manhattan +
  geom_point(color = sp_col, size = 0.4, alpha = 0.8) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(x = "position (bp)", y = "-log<sub>10</sub>(p)",
       title = if (preview_mode) paste0(opt$sp, " ", opt$trait, " GWAS: chromosome ", opt$chr_label, " (preview, all scaffolds)")
               else if (multiregion_mode) paste0(opt$sp, " ", opt$trait, " GWAS: chromosome ", opt$chr_label, " (full chromosome, all scaffolds)")
               else paste0(opt$sp, " ", opt$trait, " GWAS: ", opt$scaffold, " (full chromosome)")) +
  theme_bw() +
  theme(
    axis.title.y = element_markdown(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

ggsave(manhattan_png, plot = p_manhattan, width = 12, height = 3.5, dpi = 600, bg = "white")
message("Written: ", manhattan_png)

if (preview_mode) {
  message("Preview mode: chromosome-wide Manhattan plot written. ",
          "Pick a scaffold + start/end from this figure, then rerun with --scaffold/--start/--end for the LD heatmap.")
  quit(save = "no", status = 0)
}

# 3. LD heatmap
plink_geno <- read.plink(bed, bim, fam)
geno <- plink_geno$genotypes
map  <- plink_geno$map  # chromosome, snp.name, cM, position, allele.1, allele.2

if (multiregion_mode) {
  map$row_idx <- seq_len(nrow(map))
  map_adj <- map %>%
    inner_join(group_scaffolds %>% select(scaffold, offset), by = c("chromosome" = "scaffold")) %>%
    mutate(pos_adj = position + offset) %>%
    filter(pos_adj >= opt$start, pos_adj <= opt$end) %>%
    arrange(pos_adj)
  keep_idx <- map_adj$row_idx
  region_geno <- geno[, keep_idx]
  region_map  <- map_adj
  region_positions <- region_map$pos_adj
  ld_title_range <- paste0("chr", opt$chr_label, ":", opt$start, "-", opt$end)
} else {
  keep <- map$chromosome == opt$scaffold & map$position >= opt$start & map$position <= opt$end
  region_geno <- geno[, keep]
  region_map  <- map[keep, ]
  region_positions <- region_map$position
  ld_title_range <- paste0(opt$scaffold, ":", opt$start, "-", opt$end)
}
message("SNPs in region (PLINK data), before MAF filter: ", ncol(region_geno))

## MAF filtering: r^2
if (ncol(region_geno) > 0) {
  maf <- col.summary(region_geno)$MAF
  keep_maf <- !is.na(maf) & maf >= opt$min_maf
  message("MAF >= ", opt$min_maf, ": ", sum(keep_maf), "/", ncol(region_geno), " SNPs retained")
  region_geno      <- region_geno[, keep_maf]
  region_positions <- region_positions[keep_maf]
}
message("SNPs in region (PLINK data), after MAF filter: ", ncol(region_geno))

if (ncol(region_geno) > opt$max_ld_snps) {
  ## Evenly thin by position
  order_idx <- order(region_positions)
  thin_idx  <- order_idx[round(seq(1, length(order_idx), length.out = opt$max_ld_snps))]
  message("Thinning ", ncol(region_geno), " SNPs down to ", length(thin_idx))
  region_geno      <- region_geno[, thin_idx]
  region_positions <- region_positions[thin_idx]
}

if (ncol(region_geno) < 2) {
  warning("Fewer than 2 SNPs in region in PLINK data -- skipping LD heatmap.")
  ld_png <- NA_character_
} else {
  ld_ramp <- colorRampPalette(c(sp_col, "#FFFFEB"))(20)

  ## Compute r^2
  ld_obj <- ld(region_geno, depth = ncol(region_geno) - 1, stats = "R.squared", symmetric = TRUE)
  ld_mat <- as.matrix(ld_obj)
  off_diag_max <- max(ld_mat[upper.tri(ld_mat)], na.rm = TRUE)
  message("Max observed r^2 in this window: ", signif(off_diag_max, 3),
          " -- rescaling colour scale so this maps to the top of the ramp.")
  ld_mat_scaled <- ld_mat / off_diag_max
  diag(ld_mat_scaled) <- 1
  ld_mat_scaled[ld_mat_scaled > 1] <- 1

# putative inversion location annotations
snp_labels <- NULL
if (!is.na(opt$mark_start) || !is.na(opt$mark_end)) {
  marker_names <- colnames(ld_mat_scaled)
  if (is.null(marker_names)) {
    warning("LD matrix has no marker names -- cannot place --mark-start/--mark-end labels.")
  } else {
    marks <- c(opt$mark_start, opt$mark_end)
    marks <- marks[!is.na(marks)]
    snp_labels <- unique(sapply(marks, function(m) marker_names[which.min(abs(region_positions - m))]))
    message("Flagging nearest marker(s) to breakpoint(s) ", paste(marks, collapse = ", "), ": ",
            paste(snp_labels, collapse = ", "))
  }
}

  ld_png <- file.path(opt$outdir, paste0("LDheatmap_", tag, ".pdf"))
  cairo_pdf(ld_png, width = 16, height = 8)
  print(LDheatmap(
    ld_mat_scaled,
    genetic.distances = region_positions,
    distances = "physical",
    color = ld_ramp,
    title = paste0("LD: ", opt$sp, " ", ld_title_range, " (colour scale capped at max r\u00b2 = ", signif(off_diag_max, 3), ")"),
    add.map = !is.null(snp_labels),
    SNP.name = snp_labels,
    flip = TRUE
  ))

  dev.off()
  message("Written: ", ld_png)
}

colnames(ld_mat_scaled)[which.min(abs(region_positions - 250000))]