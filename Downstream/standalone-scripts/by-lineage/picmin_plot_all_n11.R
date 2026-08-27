#!/usr/bin/env Rscript
# plot_picmin_genome_all_species_n11.R
# PicMin genome plot for n_est=11, all species except Hana, species colours.
# Run from: gwas_pipeline/

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

trait         <- "categorical"
fdr_threshold <- 0.01

PICMIN_DIR  <- Sys.getenv("PICMIN_DIR", unset = "results/picmin")
ORTHOGROUPS <- Sys.getenv("ORTHOGROUPS_PATH",
                          unset = "/home/ss3335/rds/rds-jiggins-rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Archive/06_OrthoFinder/proteomes/OrthoFinder/Results/Orthogroups/Orthogroups.tsv")
RDS_BASE <- "/home/ss3335/rds/rds-jiggins-rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution"
DS_LONG  <- "results/driving-genes/driving_genes_ds_categorical_long.csv"

SPECIES_COLOURS <- c(
  Hera = "#3d1f5c", Hsar = "#8a6aaa",
  Hmel = "#8c3a4a", Hnum = "#c47a8a",
  Mlys = "#6eaaa0", Mpol = "#2d6b63", Mmes = "#133a36",
  Mmot = "#7a3b10", Mmen = "#c47c3a",
  Isal = "#34471f"
)

# ── Species configuration ─────────────────────────────────────────────────────
# join_scaffold : TRUE  = join on gene_clean + scaffold_acc
#                FALSE = join on gene_clean only (globally unique IDs)
# group_scaffolds: TRUE = group sub-scaffolds by chr_num (e.g. Herato0101->chr1)
# bed_to_fai_chrom: TRUE = remap BED chrom names to FAI chrom names by size rank

SPECIES_CONFIG <- list(
  Hera = list(
    of_col          = "Heliconius_erato",
    full_name       = "Heliconius erato",
    bed             = "results/snp_gene_map/genes_Hera.bed",
    fai_dir         = file.path(RDS_BASE, "Hera/03_genome"),
    chrom_pat       = "Herato",
    scaffold_re     = "(Herato[0-9]+)",
    chr_num_re      = "Herato([0-9]{2})[0-9]{2}",  # first 2 digits = chromosome
    max_chr         = 21,
    min_len         = 0,
    join_scaffold   = FALSE,
    group_scaffolds = TRUE,
    bed_to_fai_chrom = FALSE
  ),
  Hmel = list(
    of_col          = "Heliconius_melpomene",
    full_name       = "Heliconius melpomene",
    bed             = "results/snp_gene_map/genes_Hmel.bed",
    fai_dir         = file.path(RDS_BASE, "Hmel/03_genome"),
    chrom_pat       = "Hmr",
    scaffold_re     = "(Hmr[0-9]+)",
    chr_num_re      = "Hmr0*([0-9]+)",
    max_chr         = 21,
    min_len         = 9000000,
    join_scaffold   = TRUE,    # gene IDs (000160) NOT globally unique - need scaffold
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  ),
  Hnum = list(
    of_col          = "Heliconius_numata",
    full_name       = "Heliconius numata",
    bed             = "results/snp_gene_map/genes_Hnum.bed",
    fai_dir         = file.path(RDS_BASE, "Hnum/03_genome"),
    chrom_pat       = "Hnum",
    scaffold_re     = "(Hnum[0-9]+)",
    chr_num_re      = "Hnum([0-9]+)",
    max_chr         = 21,
    min_len         = 0,
    join_scaffold   = FALSE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  ),
  Hsar = list(
    of_col          = "Heliconius_sara",
    full_name       = "Heliconius sara",
    bed             = "results/snp_gene_map/genes_Hsar.bed",
    fai_dir         = file.path(RDS_BASE, "Hsar/03_genome"),
    chrom_pat       = "SUPER_",
    scaffold_re     = "(SUPER_[0-9]+)",
    chr_num_re      = "SUPER_([0-9]+)",
    max_chr         = 21,
    min_len         = 0,
    join_scaffold   = FALSE,   # g1, g2... globally unique
    group_scaffolds = FALSE,
    bed_to_fai_chrom = TRUE    # BED=HAP1_SCAFFOLD_N, FAI=SUPER_N
  ),
  Isal = list(
    of_col          = "Ithomia_salapia",
    full_name       = "Ithomia salapia",
    bed             = "results/snp_gene_map/genes_Isal.bed",
    fai_dir         = file.path(RDS_BASE, "Isal/03_genome"),
    chrom_pat       = "SUPER_",
    scaffold_re     = "(SUPER_[A-Z0-9]+)",   # matches SUPER_W, SUPER_20 etc.
    chr_num_re      = "SUPER_([0-9]+)",       # numeric only for sorting
    max_chr         = 34,
    min_len         = 0,
    join_scaffold   = TRUE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = TRUE    # BED=HAP1_SCAFFOLD_N -> remap to SUPER_N
  ),
  Mlys = list(
    of_col          = "Mechanitis_lysimnia",
    full_name       = "Mechanitis lysimnia",
    bed             = "results/snp_gene_map/genes_Mlys.bed",
    fai_dir         = file.path(RDS_BASE, "Mlys/03_genome"),
    chrom_pat       = "ENA",
    scaffold_re     = "(OZ221[0-9]+\\.[0-9]+)",
    chr_num_re      = NULL,
    max_chr         = 13,
    min_len         = 0,
    join_scaffold   = TRUE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  ),
  Mpol = list(
    of_col          = "Mechanitis_polymnia",
    full_name       = "Mechanitis polymnia",
    bed             = "results/snp_gene_map/genes_Mpol.bed",
    fai_dir         = file.path(RDS_BASE, "Mpol/03_genome"),
    chrom_pat       = "ENA",
    scaffold_re     = "(OZ240[0-9]+\\.[0-9]+)",
    chr_num_re      = NULL,
    max_chr         = 16,
    min_len         = 0,
    join_scaffold   = TRUE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  ),
  Mmes = list(
    of_col          = "Mechanitis_messenoides",
    full_name       = "Mechanitis messenoides",
    bed             = "results/snp_gene_map/genes_Mmes.bed",
    fai_dir         = file.path(RDS_BASE, "Mmes/03_genome"),
    chrom_pat       = "ENA",
    scaffold_re     = "(OY365[0-9]+\\.[0-9]+)",
    chr_num_re      = NULL,
    max_chr         = 15,
    min_len         = 0,
    join_scaffold   = TRUE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  ),
  Mmot = list(
    of_col          = "Melinaea_mothone",
    full_name       = "Melinaea mothone",
    bed             = "results/snp_gene_map/genes_Mmot.bed",
    fai_dir         = file.path(RDS_BASE, "Mmot/03_genome"),
    chrom_pat       = "ENA",
    scaffold_re     = "(OZ240[0-9]+\\.[0-9]+)",
    chr_num_re      = NULL,
    max_chr         = 14,
    min_len         = 0,
    join_scaffold   = TRUE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  ),
  Mmen = list(
    of_col          = "Melinaea_menophilus",
    full_name       = "Melinaea menophilus",
    bed             = "results/snp_gene_map/genes_Mmen.bed",
    fai_dir         = file.path(RDS_BASE, "Mmen/03_genome"),
    chrom_pat       = "ENA",
    scaffold_re     = "((?:OU911|OZ336)[0-9]+\\.[0-9]+)",
    chr_num_re      = NULL,
    max_chr         = 22,
    min_len         = 0,
    join_scaffold   = TRUE,
    group_scaffolds = FALSE,
    bed_to_fai_chrom = FALSE
  )
)

# ── Helper functions ──────────────────────────────────────────────────────────
make_extract_fn <- function(pattern) {
  function(x) {
    m <- regmatches(x, regexec(pattern, x))
    sapply(m, function(v) if (length(v) >= 2) v[2] else NA_character_)
  }
}

clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

base_theme <- theme_half_open() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13)
  )

# ── Load shared data ──────────────────────────────────────────────────────────
message("Loading PicMin results...")
picmin_file <- file.path(PICMIN_DIR, paste0("picmin_", trait, "_results.csv"))
if (!file.exists(picmin_file)) stop("Not found: ", picmin_file)
pm <- read.csv(picmin_file, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(p), n_species == 11)
message("  Orthogroups: ", nrow(pm))

message("Loading OrthoFinder table...")
og_raw <- read.delim(ORTHOGROUPS, check.names = FALSE, stringsAsFactors = FALSE)

message("Loading DS driving genes...")
ds_long <- read.csv(DS_LONG, stringsAsFactors = FALSE)

# ── Loop over species ─────────────────────────────────────────────────────────
for (sp in names(SPECIES_CONFIG)) {
  cfg <- SPECIES_CONFIG[[sp]]
  message("\n=== ", sp, " | ", cfg$full_name, " ===")
  
  out_dir <- file.path("results/picmin/by_lineage", sp)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  sp_colour    <- SPECIES_COLOURS[sp]
  extract_scaf <- make_extract_fn(cfg$scaffold_re)
  
  # ── 1. OrthoFinder gene list ────────────────────────────────────────────────
  og_long <- og_raw[, c("Orthogroup", cfg$of_col)] %>%
    setNames(c("Orthogroup", "gene_list")) %>%
    dplyr::filter(!is.na(gene_list), gene_list != "") %>%
    tidyr::separate_rows(gene_list, sep = ",\\s*") %>%
    mutate(
      gene_clean   = clean_gene_id(gene_list),
      scaffold_acc = extract_scaf(gene_list)
    ) %>%
    distinct(Orthogroup, gene_clean, .keep_all = TRUE)
  
  message("  OrthoFinder entries: ", nrow(og_long))
  
  # ── 2. Gene BED ─────────────────────────────────────────────────────────────
  if (!file.exists(cfg$bed)) {
    message("  SKIP — BED not found: ", cfg$bed); next
  }
  genes <- read.delim(cfg$bed, header = FALSE,
                      col.names = c("chrom", "start", "end", "gene_id"),
                      colClasses = c("character", "integer", "integer", "character"),
                      stringsAsFactors = FALSE) %>%
    mutate(
      gene_clean   = clean_gene_id(gene_id),
      scaffold_acc = extract_scaf(chrom)
    )
  
  message("  BED entries: ", nrow(genes))
  
  # ── 3. Load FAI early (needed for bed_to_fai_chrom and scaffold grouping) ───
  fai_files <- list.files(cfg$fai_dir,
                          pattern = "\\.fa\\.fai$|\\.fasta\\.fai$|\\.fna\\.fai$",
                          full.names = TRUE)
  if (length(fai_files) == 0) {
    message("  SKIP — no FAI in: ", cfg$fai_dir); next
  }
  fai <- read.delim(fai_files[1], header = FALSE)
  colnames(fai)[1:2] <- c("scaffold", "length")
  
  main <- fai %>%
    dplyr::filter(grepl(cfg$chrom_pat, scaffold, fixed = FALSE),
                  length >= cfg$min_len)
  
  if (nrow(main) == 0) {
    message("  SKIP — no scaffolds matching: ", cfg$chrom_pat); next
  }
  
  # ── 4. Remap BED chroms to FAI chroms if needed ──────────────────────────────
  if (isTRUE(cfg$bed_to_fai_chrom)) {
    fai_ordered <- main %>% arrange(desc(length)) %>%
      mutate(rank = seq_len(n()))
    
    bed_chrom_sizes <- genes %>%
      group_by(chrom) %>%
      summarise(total_len = max(end), .groups = "drop") %>%
      arrange(desc(total_len)) %>%
      mutate(rank = seq_len(n()))
    
    chrom_map <- bed_chrom_sizes %>%
      left_join(fai_ordered %>% dplyr::select(rank, fai_chrom = scaffold),
                by = "rank")
    
    genes <- genes %>%
      left_join(chrom_map %>% dplyr::select(chrom, fai_chrom), by = "chrom") %>%
      mutate(
        scaffold_acc = ifelse(!is.na(fai_chrom), extract_scaf(fai_chrom), scaffold_acc),
        chrom        = ifelse(!is.na(fai_chrom), fai_chrom, chrom)
      ) %>%
      dplyr::select(-fai_chrom)
    
    message("  BED chroms remapped to FAI, sample: ",
            paste(head(unique(genes$chrom), 3), collapse = ", "))
  }
  
  # ── 5. og_pos — join OrthoFinder genes to BED positions ─────────────────────
  if (isTRUE(cfg$join_scaffold)) {
    og_pos_all <- og_long %>%
      inner_join(genes, by = c("gene_clean", "scaffold_acc"),
                 relationship = "many-to-many") %>%
      group_by(Orthogroup) %>% slice_head(n = 1) %>% ungroup() %>%
      mutate(pos = (start + end) / 2) %>%
      dplyr::select(Orthogroup, chrom, pos)
  } else {
    og_pos_all <- og_long %>%
      inner_join(genes %>% dplyr::select(gene_clean, chrom, start, end),
                 by = "gene_clean", relationship = "many-to-many") %>%
      group_by(Orthogroup) %>% slice_head(n = 1) %>% ungroup() %>%
      mutate(pos = (start + end) / 2) %>%
      dplyr::select(Orthogroup, chrom, pos)
  }
  
  message("  og_pos_all: ", nrow(og_pos_all))
  
  # DS driving gene override for significant OGs
  ds_sp <- ds_long %>%
    dplyr::filter(species == sp) %>%
    dplyr::select(Orthogroup, driving_gene_raw) %>%
    mutate(
      gene_clean   = clean_gene_id(driving_gene_raw),
      scaffold_acc = extract_scaf(driving_gene_raw)
    )
  
  if (nrow(ds_sp) > 0) {
    if (isTRUE(cfg$join_scaffold)) {
      ds_pos <- ds_sp %>%
        inner_join(genes, by = c("gene_clean", "scaffold_acc"),
                   relationship = "many-to-many") %>%
        group_by(Orthogroup) %>% slice_head(n = 1) %>% ungroup() %>%
        mutate(pos = (start + end) / 2) %>%
        dplyr::select(Orthogroup, chrom, pos)
    } else {
      ds_pos <- ds_sp %>%
        inner_join(genes %>% dplyr::select(gene_clean, chrom, start, end),
                   by = "gene_clean", relationship = "many-to-many") %>%
        group_by(Orthogroup) %>% slice_head(n = 1) %>% ungroup() %>%
        mutate(pos = (start + end) / 2) %>%
        dplyr::select(Orthogroup, chrom, pos)
    }
    og_pos <- og_pos_all %>%
      dplyr::filter(!Orthogroup %in% ds_pos$Orthogroup) %>%
      bind_rows(ds_pos)
    message("  DS positions: ", nrow(ds_pos))
  } else {
    og_pos <- og_pos_all
  }
  
  message("  Total OGs with positions: ", nrow(og_pos))
  
  # ── 6. Build chromosome layout ───────────────────────────────────────────────
  if (!is.null(cfg$chr_num_re)) {
    m_chr <- regmatches(main$scaffold, regexec(cfg$chr_num_re, main$scaffold))
    num_extract <- suppressWarnings(
      as.integer(sapply(m_chr, function(v) if (length(v) >= 2) v[2] else NA))
    )
    main$chr_num <- num_extract
  } else {
    main$chr_num <- NA_integer_
  }
  
  if (isTRUE(cfg$group_scaffolds) && sum(!is.na(main$chr_num)) > 0) {
    # Group sub-scaffolds by chromosome number (e.g. Herato0101, Herato0102 -> chr 1)
    chr_ranks <- data.frame(
      chr_num = sort(unique(na.omit(main$chr_num)))
    ) %>% mutate(CHR = seq_len(n()))
    
    main <- main %>%
      left_join(chr_ranks, by = "chr_num") %>%
      dplyr::filter(!is.na(CHR), CHR <= cfg$max_chr) %>%
      arrange(CHR, desc(length))
    
    # Global offset across all sub-scaffolds
    main$global_offset <- cumsum(lag(main$length, default = 0))
    
    scaffold_to_chr <- main %>%
      dplyr::select(chrom = scaffold, CHR, length, global_offset)
    
    chr_layout <- scaffold_to_chr %>%
      group_by(CHR) %>%
      summarise(
        xmin   = min(global_offset),
        xmax   = max(global_offset + length),
        length = sum(length),
        offset = min(global_offset),
        .groups = "drop"
      )
    
  } else {
    # One scaffold per chromosome — sort and assign CHR
    if (sum(!is.na(main$chr_num)) > nrow(main) * 0.5) {
      main <- main %>% arrange(chr_num)
    } else {
      main <- main %>% arrange(desc(length))
    }
    main$CHR <- seq_len(nrow(main))
    main <- main %>% dplyr::filter(CHR <= cfg$max_chr)
    main$global_offset <- cumsum(lag(main$length, default = 0))
    
    scaffold_to_chr <- main %>%
      dplyr::select(chrom = scaffold, CHR, length, global_offset)
    
    chr_layout <- scaffold_to_chr %>%
      mutate(xmin = global_offset, xmax = global_offset + length,
             offset = global_offset)
  }
  
  message("  Chromosomes used: ", n_distinct(scaffold_to_chr$CHR))
  
  # ── 7. Build og_genome ───────────────────────────────────────────────────────
  unplaced_num <- max(scaffold_to_chr$CHR) + 1L
  
  og_genome <- og_pos %>%
    left_join(scaffold_to_chr %>% dplyr::select(chrom, CHR, global_offset),
              by = "chrom") %>%
    mutate(
      CHR           = ifelse(is.na(CHR), unplaced_num, CHR),
      global_offset = ifelse(is.na(global_offset),
                             max(scaffold_to_chr$global_offset +
                                   scaffold_to_chr$length), global_offset),
      x_pos = pos + global_offset
    ) %>%
    left_join(pm %>% dplyr::select(Orthogroup, p, q, n_est), by = "Orthogroup") %>%
    dplyr::filter(!is.na(p), CHR != unplaced_num) %>%
    mutate(neg_log_p = -log10(p))
  
  message("  Orthogroups plotted: ", nrow(og_genome))
  
  # ── 8. Labels and bands ──────────────────────────────────────────────────────
  chr_midpoints <- chr_layout %>%
    mutate(mid = (xmin + xmax) / 2, label = as.character(CHR))
  
  chr_bands <- chr_layout %>%
    mutate(fill_col = ifelse(CHR %% 2 == 0, "grey92", "white"))
  
  band_annotations <- lapply(seq_len(nrow(chr_bands)), function(i) {
    annotate("rect",
             xmin = chr_bands$xmin[i], xmax = chr_bands$xmax[i],
             ymin = -Inf, ymax = Inf,
             fill = chr_bands$fill_col[i], alpha = 0.6)
  })
  
  # ── 9. Plot n_est=11 ─────────────────────────────────────────────────────────
  highlight  <- og_genome %>% dplyr::filter(n_est == 11, q < fdr_threshold)
  background <- og_genome %>% dplyr::filter(!(n_est == 11 & q < fdr_threshold))
  message("  Highlighted OGs (n_est=11): ", nrow(highlight))
  
  p <- ggplot() +
    band_annotations +
    geom_point(data = background,
               aes(x = x_pos, y = neg_log_p),
               colour = "grey60", alpha = 0.2, size = 1.6) +
    geom_point(data = highlight,
               aes(x = x_pos, y = neg_log_p),
               fill = sp_colour, colour = "grey20",
               shape = 21, size = 2.8, alpha = 0.9, stroke = 0.3) +
    geom_hline(yintercept = -log10(fdr_threshold),
               lty = 2, colour = "grey30", linewidth = 0.5) +
    scale_x_continuous(
      breaks = chr_midpoints$mid,
      labels = chr_midpoints$label,
      expand = c(0.01, 0.01)
    ) +
    labs(
      x     = "Chromosome",
      y     = expression(-log[10]*"(PicMin p-value)"),
      title = paste0("PicMin — Categorical | ", cfg$full_name,
                     " | n_est = 11 (n=", nrow(highlight), " orthogroups)")
    ) +
    base_theme
  
  out_png <- file.path(out_dir,
                       paste0("picmin_", trait, "_", sp, "_genome_n11.png"))
  out_pdf <- file.path(out_dir,
                       paste0("picmin_", trait, "_", sp, "_genome_n11.pdf"))
  
  ggsave(out_png, p, width = 16, height = 6, dpi = 250, bg = "white")
  ggsave(out_pdf, p, width = 16, height = 6, bg = "white", device = cairo_pdf)
  message("  Saved: ", basename(out_png))
}

message("\nDone.")