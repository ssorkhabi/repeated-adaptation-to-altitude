#!/usr/bin/env Rscript
# manhattan.R (ggplot2)
#
# Genome-wide Manhattan plot

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(ggplot2)
  library(ggtext)
})

## 0. Species colours
this_script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(this_script_dir) == 0 || this_script_dir == "") this_script_dir <- "."
source(file.path(this_script_dir, "species_config.R"))

## 1. Args
option_list <- list(
  make_option("--assoc",         type = "character", help = "GEMMA assoc file"),
  make_option("--genome-dir",    type = "character", dest = "genome_dir", help = "Directory containing the .fai"),
  make_option("--out",           type = "character", help = "Output PNG path"),
  make_option("--sp",            type = "character", help = "Species code (must be a key in species_colors)"),
  make_option("--trait",         type = "character", default = "altitude", help = "Trait label for plot title [default %default]"),
  make_option("--chrom-pattern", type = "character", dest = "chrom_pat", help = "Regex to select main scaffolds, e.g. '^Herato' or '^SUPER_' or '^ENA'"),
  make_option("--max-chr",       type = "integer",   default = NA_integer_, dest = "max_chr", help = "Restrict to first N autosomes (sex chromosomes always kept)"),
  make_option("--sex-chr-override", type = "character", default = NA_character_, dest = "sex_override",
              help = "Comma-separated CHR:label overrides for sex chromosomes not catchable by name pattern, e.g. '21:Z'"),
  make_option("--sex-pattern",   type = "character", default = "(_Z|_W)$", dest = "sex_pattern", help = "Regex (applied to scaffold name) to auto-detect sex scaffolds [default %default]"),
  make_option("--regions",       type = "character", default = NA_character_, help = "TSV of candidate regions to highlight (scaffold, start, end, label)"),
  make_option("--pval-col",      type = "character", default = "p_wald", dest = "pval_col", help = "P-value column [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (req in c("assoc", "genome_dir", "out", "sp", "chrom_pat")) {
  if (is.null(opt[[req]])) stop("--", gsub("_", "-", req), " is required")
}
if (!opt$sp %in% names(species_colors)) {
  stop("--sp must be one of: ", paste(names(species_colors), collapse = ", "))
}
sp_col <- species_colors[[opt$sp]]

## 2. Read and thin GWAS results
gwas <- read.delim(opt$assoc, sep = "\t", header = TRUE)
gwas <- gwas[!grepl("^##", gwas$chr), ]
gwas <- gwas[!is.na(gwas[[opt$pval_col]]), ]
gwas$P <- gwas[[opt$pval_col]]

n_total_tests <- nrow(gwas)
bonferroni_threshold <- 0.05 / n_total_tests
message("Total SNPs tested: ", n_total_tests, " | Bonferroni threshold: ", signif(bonferroni_threshold, 3))

sig_mask    <- gwas$P < 0.01
n_nonsig    <- sum(!sig_mask)
keep_nonsig <- sample(which(!sig_mask), size = min(n_nonsig, ceiling(n_nonsig * 0.1)))
gwas        <- gwas[c(which(sig_mask), keep_nonsig), ]
gwas        <- gwas[order(gwas$chr, gwas$ps), ]
message("SNPs after thinning: ", nrow(gwas))

## 3. Read FASTA index and build scaffold -> CHR-label mapping
fai_files <- list.files(opt$genome_dir,
  pattern = "\\.fa\\.fai$|\\.fasta\\.fai$|\\.fna\\.fai$", full.names = TRUE)
if (length(fai_files) == 0) stop("No .fai file found in ", opt$genome_dir)

fai <- read.delim(fai_files[1], header = FALSE)
colnames(fai)[1:2] <- c("scaffold", "length")

main <- fai %>% filter(grepl(opt$chrom_pat, scaffold))
if (nrow(main) == 0) stop("No scaffolds matching pattern '", opt$chrom_pat, "' found in FAI")
message("Scaffolds matching '", opt$chrom_pat, "': ", nrow(main))

## 3a. Auto-detect sex scaffolds by name (catches Isal's SUPER_Z/SUPER_W etc.)
main$is_sex_auto <- grepl(opt$sex_pattern, main$scaffold, ignore.case = TRUE)

## 3b. Decide numbering mode
non_sex_main <- main %>% filter(!is_sex_auto)

suffix_match <- regmatches(non_sex_main$scaffold, regexec("^(.*?)([0-9]+)([A-Za-z]*)$", non_sex_main$scaffold))
non_sex_main$chr_num_str <- sapply(suffix_match, function(x) if (length(x) == 4) x[3] else NA_character_)
non_sex_main$sub_letter  <- sapply(suffix_match, function(x) if (length(x) == 4) x[4] else NA_character_)
## "Clean" means the full trailing digit run is exactly 2 digits (Hmr01, Hnum05a).
## Using [0-9]+ (not [0-9]{2}) to capture the WHOLE run first, then checking its
## length here, matters: a fixed {2} capture lets the lazy prefix absorb extra
## leading digits of a longer run (e.g. Hera's Herato0201 would wrongly match
## by letting prefix eat "Herato02" and capture just "01").
has_clean_num <- !is.na(non_sex_main$chr_num_str) & nchar(non_sex_main$chr_num_str) == 2
frac_clean    <- if (nrow(non_sex_main) > 0) mean(has_clean_num) else 0
n_clean_chr   <- n_distinct(non_sex_main$chr_num_str[has_clean_num])
plausible_cap <- if (!is.na(opt$max_chr)) opt$max_chr * 2 else 40
clean_mode_ok <- nrow(non_sex_main) > 0 && frac_clean > 0.8 && n_clean_chr <= plausible_cap
if (nrow(non_sex_main) > 0 && frac_clean > 0.8 && !clean_mode_ok) {
  message("Clean-numbering mode rejected: ", n_clean_chr, " distinct chromosome numbers implied, ",
          "exceeds plausibility cap of ", plausible_cap, " -- --chrom-pattern is probably sweeping in unplaced fragments. ",
          "Falling through to grouping/length-rank mode instead.")
}

re_group   <- paste0("^", opt$chrom_pat, "(\\d{2})(\\d+)$")
n_matching <- sum(grepl(re_group, non_sex_main$scaffold))
frac_match <- ifelse(nrow(non_sex_main) > 0, n_matching / nrow(non_sex_main), 0)
message("Scaffolds matching grouping regex: ", n_matching, "/", nrow(non_sex_main),
        " (", round(frac_match * 100), "%)")

if (clean_mode_ok) {
  message("Clean-numbering mode: using each scaffold's own embedded chromosome number ",
          "(handles both one-scaffold-per-chromosome, e.g. Hmr01, and lettered sub-scaffolds sharing a number, e.g. Hnum05a/b/c)")
  autosomes <- non_sex_main %>%
    filter(has_clean_num) %>%
    mutate(CHR_num = as.integer(chr_num_str), CHR_label = as.character(CHR_num), chr_group = NA_character_) %>%
    arrange(CHR_num, sub_letter, scaffold)
} else if (n_matching > 5 && frac_match > 0.5) {
  message("Grouping mode: extracting embedded chromosome number")
  chr_extracted <- regmatches(non_sex_main$scaffold, regexec(re_group, non_sex_main$scaffold))
  non_sex_main$chr_group <- sapply(chr_extracted, function(x) if (length(x) == 3) x[2] else NA_character_)
  grouped <- non_sex_main %>%
    filter(!is_sex_auto, !is.na(chr_group)) %>%
    arrange(as.integer(chr_group), scaffold)
  chr_nums <- sort(unique(as.integer(grouped$chr_group)))

  if (length(chr_nums) >= 5 && length(chr_nums) < nrow(non_sex_main)) {
    grouped$CHR_num   <- match(as.integer(grouped$chr_group), chr_nums)
    grouped$CHR_label <- as.character(grouped$CHR_num)
    autosomes <- grouped
  } else {
    message("  Grouping produced too few/many groups -- falling back to per-scaffold mode")
    autosomes <- non_sex_main %>%
      arrange(desc(length)) %>%
      mutate(CHR_num = row_number(), CHR_label = as.character(CHR_num), chr_group = NA_character_)
  }
} else {
  message("Per-scaffold mode (sorted by size)")
  autosomes <- non_sex_main %>%
    arrange(desc(length)) %>%
    mutate(CHR_num = row_number(), CHR_label = as.character(CHR_num), chr_group = NA_character_)
}

## 3c. Sex scaffolds
sex_main <- main %>% filter(is_sex_auto) %>% arrange(desc(length))
if (nrow(sex_main) > 0) {
  sex_main$CHR_label <- toupper(sub(".*_([A-Za-z]+)$", "\\1", sex_main$scaffold))
  sex_main$chr_group <- NA_character_
  sex_main$CHR_num   <- max(autosomes$CHR_num, 0) + seq_len(nrow(sex_main))
} else {
  sex_main$CHR_label <- character(0)
  sex_main$chr_group <- character(0)
  sex_main$CHR_num   <- integer(0)
}

main2 <- bind_rows(
  autosomes %>% select(scaffold, length, chr_group, CHR_num, CHR_label),
  sex_main  %>% select(scaffold, length, chr_group, CHR_num, CHR_label)
)

## 3d. Manual override, e.g. --sex-chr-override "21:Z" for species where the
##     sex chromosome isn't name-distinguishable (Heliconius chr21 = Z convention)
override_labels <- character(0)
if (!is.na(opt$sex_override)) {
  pairs <- strsplit(strsplit(opt$sex_override, ",")[[1]], ":")
  for (p in pairs) {
    chr_to_relabel <- trimws(p[1])
    new_label      <- trimws(p[2])
    n_affected <- sum(main2$CHR_label == chr_to_relabel)
    message("Overriding CHR '", chr_to_relabel, "' -> '", new_label, "' (", n_affected, " scaffolds)")
    main2$CHR_label[main2$CHR_label == chr_to_relabel] <- new_label
    override_labels <- c(override_labels, new_label)
  }
}

## 3e. Optionally restrict autosome count (sex chromosomes always kept --
##     both auto-detected ones AND manually-overridden ones, e.g. Hera's
##     chr21->Z, which wouldn't otherwise appear in sex_main at all)
if (!is.na(opt$max_chr)) {
  keep_labels <- c(as.character(seq_len(opt$max_chr)), unique(sex_main$CHR_label), override_labels)
  n_before <- n_distinct(main2$CHR_label)
  main2 <- main2 %>% filter(CHR_label %in% keep_labels)
  message("Restricting to first ", opt$max_chr, " autosomes + sex chromosomes (discarded ",
          n_before - n_distinct(main2$CHR_label), " groups)")
}

## 3f. Cumulative offset WITHIN each CHR_label (fixes overlapping scaffolds
##     in grouping mode; a no-op for per-scaffold mode where each CHR_label
##     has exactly one scaffold)
main2 <- main2 %>%
  arrange(CHR_num, scaffold) %>%
  group_by(CHR_label) %>%
  mutate(offset = lag(cumsum(length), default = 0)) %>%
  ungroup()

chr_lengths <- main2 %>%
  group_by(CHR_label, CHR_num) %>%
  summarise(chr_length = sum(length), .groups = "drop") %>%
  arrange(CHR_num)

## Save the lookup table -- this is the artifact that removes all future
## ambiguity about what a plotted chromosome number/label actually is.
lookup_out <- sub("\\.png$", "_chr_lookup.csv", opt$out)
dir.create(dirname(lookup_out), showWarnings = FALSE, recursive = TRUE)
write.csv(main2 %>% select(scaffold, length, CHR_label, CHR_num, offset),
          lookup_out, row.names = FALSE)
message("Written chromosome lookup table: ", lookup_out)

## 4. Join CHR_label + adjusted BP into GWAS data
gwas <- gwas %>%
  inner_join(main2 %>% select(scaffold, CHR_label, CHR_num, offset),
             by = c("chr" = "scaffold")) %>%
  mutate(BP_adj = ps + offset) %>%
  filter(!is.na(P), P > 0)

message("SNPs retained on placed chromosomes: ", nrow(gwas))
message("Chromosomes present: ", paste(sort(unique(gwas$CHR_label)), collapse = ", "))

## 5. Genome-wide cumulative position (bp_cum)
dup_labels <- chr_lengths %>% count(CHR_label) %>% filter(n > 1) %>% pull(CHR_label)
if (length(dup_labels) > 0) {
  stop("CHR_label(s) [", paste(dup_labels, collapse = ", "), "] map to more than one distinct group. ",
       "This usually means --sex-chr-override relabeled a scaffold to a name that's already ",
       "auto-detected by --sex-pattern (e.g. the species already has a named SUPER_Z/SUPER_W). ",
       "Check the chr_lookup.csv and drop the override if the sex chromosome is already name-detectable.")
}

chr_lengths <- chr_lengths %>%
  arrange(CHR_num) %>%
  mutate(genome_offset = lag(cumsum(chr_length), default = 0),
         mid = genome_offset + chr_length / 2,
         plot_rank = row_number())

gwas <- gwas %>%
  left_join(chr_lengths %>% select(CHR_label, genome_offset, plot_rank), by = "CHR_label") %>%
  mutate(BP_genome = BP_adj + genome_offset)

## Alternating colours by chromosome (odd/even), same visual role as the
## classic two-tone qqman look, built from the species' own colour.
light_col <- colorRampPalette(c(sp_col, "white"))(5)[2]
dark_col  <- colorRampPalette(c(sp_col, "black"))(5)[2]
gwas$point_col <- ifelse(gwas$plot_rank %% 2 == 0, light_col, dark_col)

## 6. Optional candidate-region highlights (raw scaffold coords -> adjusted)
highlight_df <- NULL
if (!is.na(opt$regions) && file.exists(opt$regions)) {
  regions <- read.delim(opt$regions, header = TRUE)
  regions <- regions %>%
    inner_join(main2 %>% select(scaffold, CHR_label, offset), by = "scaffold") %>%
    inner_join(chr_lengths %>% select(CHR_label, genome_offset), by = "CHR_label") %>%
    mutate(start_adj = start + offset + genome_offset, end_adj = end + offset + genome_offset)
  if (nrow(regions) > 0) highlight_df <- regions
  message("Loaded ", nrow(regions), " candidate region(s) to highlight")
}

## 7. Plot
ymax <- ceiling(max(-log10(gwas$P), na.rm = TRUE)) + 1

## Alternating white/grey background bands
band_df <- chr_lengths %>%
  mutate(xmin = genome_offset, xmax = genome_offset + chr_length,
         band_col = ifelse(plot_rank %% 2 == 0, "grey92", "white"))

p <- ggplot(gwas, aes(x = BP_genome, y = -log10(P))) +
  geom_rect(data = band_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_col),
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_point(aes(color = point_col), size = 0.4, alpha = 0.8) +
  scale_color_identity() +
  scale_x_continuous(breaks = chr_lengths$mid, labels = chr_lengths$CHR_label,
                      expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0.02, 0.02))) +
  geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  labs(x = NULL, y = NULL,
       title = paste0(opt$sp, " ", opt$trait, " GWAS")) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(size = 10),
    axis.text.x = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

if (!is.null(highlight_df)) {
  p <- p + geom_rect(data = highlight_df,
                      aes(xmin = start_adj, xmax = end_adj, ymin = 0, ymax = ymax),
                      inherit.aes = FALSE, fill = sp_col, alpha = 0.18)
}

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
plot_width <- max(12, nrow(chr_lengths) * 0.4)
ggsave(opt$out, plot = p, width = plot_width, height = 5, dpi = 600, bg = "white")
message("Written: ", opt$out)