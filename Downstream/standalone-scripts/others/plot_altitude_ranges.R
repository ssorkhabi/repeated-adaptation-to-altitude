#!/usr/bin/env Rscript
# plot_altitude_ranges.R

library(tidyverse)
library(cowplot)

# Config
PHENO_DIR  <- "phenotypes"
OUT_DIR    <- "results/plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SPECIES <- c("Mlys", "Mpol", "Mmes", "Hana", "Isal", "Mmen", "Mmot",
             "Hnum", "Hmel", "Hsar", "Hera")

SPECIES_LABELS <- c(
  Hera = "He. erato",
  Hmel = "He. melpomene",
  Hnum = "He. numata",
  Hsar = "He. sara",
  Hana = "Hy. anastasia",
  Isal = "I. salapia",
  Mlys = "Mec. lysimnia",
  Mpol = "Mec. polymnia",
  Mmes = "Mec. messenoides",
  Mmot = "Mel. mothone",
  Mmen = "Mel. menophilus"
)

# Colour palette
SP_COLS <- c(
  # Heliconius erato + sara — purples
  Hera    = "#3d1f5c",
  Hsar    = "#8a6aaa",

  # Heliconius melpomene + numata — pinks
  Hmel   = "#8c3a4a",
  Hnum   = "#c47a8a",

  # Mechanitis — blues
  Mlys   = "#6eaaa0",
  Mpol   = "#2d6b63",
  Mmes   = "#133a36",

  # Melinaea — oranges
  Mmot   = "#7a3b10",
  Mmen   = "#c47c3a",

  # Ithomia + Hypothyris — greens
  Isal   = "#34471f",
  Hana   = "#5a7a35"
)

# Read phenotype files
pheno_list <- list()

for (sp in SPECIES) {
  f <- file.path(PHENO_DIR, paste0(sp, "_phenotypes.txt"))
  if (!file.exists(f)) {
    message("WARNING: missing phenotype file for ", sp, " — skipping")
    next
  }
  d <- read.table(f, header = FALSE,
                  col.names = c("FID", "IID", "altitude"),
                  stringsAsFactors = FALSE)
  d$species     <- sp
  d$species_lab <- SPECIES_LABELS[sp]
  pheno_list[[sp]] <- d
}

pheno <- bind_rows(pheno_list) %>%
  filter(!is.na(altitude)) %>%
  mutate(
    species     = factor(species,     levels = SPECIES),
    species_lab = factor(species_lab, levels = SPECIES_LABELS[SPECIES])
  )

# Compute IQR thresholds (Q25 and Q75 per species — matches binary GWAS)
thresholds <- pheno %>%
  group_by(species, species_lab) %>%
  summarise(
    q25    = quantile(altitude, 0.25),
    q75    = quantile(altitude, 0.75),
    median = median(altitude),
    min    = min(altitude),
    max    = max(altitude),
    n      = n(),
    .groups = "drop"
  )

message("Species summary:")
print(thresholds)

# Shared theme
base_theme <- theme_half_open() +
  background_grid(major = "y") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "italic", size = 18),
    legend.position = "none"
  )

# all species side by side
p1 <- ggplot(pheno, aes(x = species_lab, y = altitude, colour = species)) +
  
  # full range line
  geom_linerange(
    data = thresholds,
    aes(x = species_lab, ymin = min, ymax = max, colour = species),
    linewidth = 0.6, alpha = 0.5, inherit.aes = FALSE
  ) +
  
  # IQR box (Q25–Q75)
  geom_rect(
    data = thresholds,
    aes(xmin = as.numeric(species_lab) - 0.25,
        xmax = as.numeric(species_lab) + 0.25,
        ymin = q25, ymax = q75,
        fill = species),
    colour = NA, alpha = 0.3, inherit.aes = FALSE
  ) +
  
  # median line
  geom_segment(
    data = thresholds,
    aes(x    = as.numeric(species_lab) - 0.25,
        xend = as.numeric(species_lab) + 0.25,
        y    = median, yend = median,
        colour = species),
    linewidth = 1, inherit.aes = FALSE
  ) +
  
  # raw points (jittered)
  geom_jitter(
    aes(fill = species),
    shape = 21, size = 1.8, alpha = 0.6,
    width = 0.18, stroke = 0.2, colour = "white"
  ) +

  scale_colour_manual(values = SP_COLS) +
  scale_fill_manual(values = SP_COLS) +
  scale_x_discrete(labels = SPECIES_LABELS[SPECIES]) +
  labs(
    x     = NULL,
    y     = "Altitude (m)",
    title = NULL,
    subtitle = NULL
  ) +
  base_theme

p1

ggsave(
  file.path(OUT_DIR, "altitude_ranges_combined.pdf"),
  p1, width = 20, height = 7, bg = "white", limitsize = FALSE,
  device = cairo_pdf
)

message("Saved: ", file.path(OUT_DIR, "altitude_ranges_combined.png"))

write.csv(pheno,file.path(OUT_DIR, "altitude_combined.csv"), row.names = FALSE)