#!/usr/bin/env Rscript
# flybase_annotation.R
#
# For each PicMin-significant orthogroup, retrieves functional annotation
# for the top Drosophila gene hit using the NCBI Drosophila gene_info file.
# No API calls needed — all data comes from the predownloaded gene_info file.
#
# Download the gene_info file first (on a login node):
#   wget "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Drosophila_melanogaster.gene_info.gz" \
#       -O results/flybase_annotation/dmel_gene_info.tsv.gz
#   gunzip results/flybase_annotation/dmel_gene_info.tsv.gz
#
# Usage: Rscript flybase_annotation.R <og_blast_summary> <gene_info> <out_dir>

suppressPackageStartupMessages({
  library(tidyverse)
})

args          <- commandArgs(trailingOnly = TRUE)
blast_summary <- args[1]
gene_info     <- args[2]
out_dir       <- if (length(args) >= 3) args[3] else "results/flybase_annotation"
# tax_id         <- 7227  # Drosophila melanogaster

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load OG blast summary
message("Loading OG blast summary...")
og_summary <- read.csv(blast_summary, stringsAsFactors = FALSE) %>%
  filter(!is.na(top_dmel_symbol), top_dmel_symbol != "unknown")

message("Orthogroups to annotate: ", nrow(og_summary))
unique_symbols <- unique(og_summary$top_dmel_symbol)
message("Unique Drosophila symbols: ", length(unique_symbols))

# 2. Load NCBI gene_info file
message("Loading NCBI Drosophila gene_info...")
if (!file.exists(gene_info)) stop("gene_info file not found: ", gene_info)

gi <- read.delim(gene_info, comment.char = "#",
                 stringsAsFactors = FALSE,
                 quote = "") %>%
  filter(tax_id == 7227) %>%
  select(
    symbol      = Symbol,
    gene_name   = Full_name_from_nomenclature_authority,
    description = description,
    synonyms    = Synonyms,
    dbXrefs     = dbXrefs
  ) %>%
  mutate(
    # Extract FlyBase FBgn ID from dbXrefs column
    fbgn_id = sub(".*FLYBASE:(FBgn[0-9]+).*", "\\1", dbXrefs),
    fbgn_id = ifelse(grepl("^FBgn", fbgn_id), fbgn_id, NA_character_),
    flybase_url = ifelse(!is.na(fbgn_id),
                         paste0("https://flybase.org/gene/", fbgn_id),
                         NA_character_)
  )

message("  Gene entries loaded: ", nrow(gi))

# 3. Match symbols to gene_info entries
# Try direct symbol match first
symbol_map <- data.frame(query_symbol = unique_symbols, stringsAsFactors = FALSE) %>%
  left_join(gi %>% select(symbol, gene_name, description, fbgn_id,
                           flybase_url, synonyms),
            by = c("query_symbol" = "symbol"))

# For unmatched: try matching against synonyms column
unmatched <- symbol_map %>% filter(is.na(fbgn_id)) %>% pull(query_symbol)
if (length(unmatched) > 0) {
  message("  Trying synonym match for ", length(unmatched), " unmatched symbols...")
  for (sym in unmatched) {
    syn_match <- gi %>%
      filter(grepl(paste0("(^|\\|)", sym, "(\\||$)"), synonyms)) %>%
      slice_head(n = 1)
    if (nrow(syn_match) > 0) {
      idx <- which(symbol_map$query_symbol == sym)
      symbol_map$gene_name[idx]    <- syn_match$gene_name
      symbol_map$description[idx]  <- syn_match$description
      symbol_map$fbgn_id[idx]      <- syn_match$fbgn_id
      symbol_map$flybase_url[idx]  <- syn_match$flybase_url
    }
  }
}

message("Symbols matched: ", sum(!is.na(symbol_map$fbgn_id)),
        " / ", nrow(symbol_map))
still_unmatched <- symbol_map$query_symbol[is.na(symbol_map$fbgn_id)]
if (length(still_unmatched) > 0) {
  message("Still unmatched: ", paste(still_unmatched, collapse = ", "))
}

# 4. Join with OG summary
final_table <- og_summary %>%
  left_join(symbol_map %>% select(-synonyms),
            by = c("top_dmel_symbol" = "query_symbol")) %>%
  select(
    Orthogroup, q, p, n_est,
    top_dmel_symbol, pct_consensus, n_species_with_hit,
    fbgn_id, flybase_url, gene_name, description, top5_symbols
  ) %>%
  arrange(q)

write.csv(final_table,
          file.path(out_dir, "og_flybase_annotation.csv"),
          row.names = FALSE)
message("Saved: og_flybase_annotation.csv")

# 5. Print summary
message("\nAnnotation of significant orthogroups!!!\n")
print(
  final_table %>%
    select(Orthogroup, q, n_est, top_dmel_symbol, gene_name, description) %>%
    mutate(description = substr(description, 1, 80)) %>%
    as.data.frame(),
  row.names = FALSE
)

message("\nDone. Full table: ",
        file.path(out_dir, "og_flybase_annotation.csv"))
