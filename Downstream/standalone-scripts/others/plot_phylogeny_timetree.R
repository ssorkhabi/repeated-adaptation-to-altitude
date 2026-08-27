#!/usr/bin/env Rscript
# plot_phylogeny_timetree.R
#
# Time-scales the OrthoFinder species tree using ape::chronos() with
# calibration points from TimeTree (timetree.org).
# Each species tip is coloured individually
#
# Node ages shown for high-confidence nodes only (from TimeTree):
#   Root (Heliconiini/Ithomiini): 77 Ma (69.4-111.5)
#   Heliconiini crown:            11.55 Ma (6.39-13.77)
#   Ithomiini crown:              24.14 Ma (22.88-25.40)


library(ape)
library(devtools)
devtools::install_github("GuangchuangYu/ggtree")
library(ggimage)
library(tidyverse)
library(cowplot)

args      <- commandArgs(trailingOnly = TRUE)
tree_file <- "/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Archive/06_OrthoFinder/proteomes/OrthoFinder/Results/Species_Tree/SpeciesTree_rooted.txt"
out_dir   <- "results/plots/phylogeny"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Species and branch colours
SP_COLOURS <- c(
  # Heliconius erato + sara — purples
  Heliconius_erato       = "#3d1f5c",
  Heliconius_sara        = "#8a6aaa",

  # Heliconius melpomene + numata — pinks
  Heliconius_melpomene   = "#8c3a4a",
  Heliconius_numata      = "#c47a8a",

  # Mechanitis — blues
  Mechanitis_lysimnia    = "#6eaaa0",
  Mechanitis_polymnia    = "#2d6b63",
  Mechanitis_messenoides = "#133a36",

  # Melinaea — oranges
  Melinaea_mothone       = "#7a3b10",
  Melinaea_menophilus    = "#c47c3a",

  # Ithomia + Hypothyris — greens
  Ithomia_salapia        = "#34471f",
  Hypothyris_anastasia   = "#5a7a35"
)

# 2. Read and time-scale tree
message("Reading tree...")
tree <- read.tree(tree_file)

message("Time-scaling with ape::chronos() using TimeTree calibrations...")
cal <- makeChronosCalib(tree,
  node = c(
    MRCA(tree, c("Heliconius_erato",  "Melinaea_mothone")),    # root
    MRCA(tree, c("Heliconius_erato",  "Heliconius_numata")),   # Heliconiini
    MRCA(tree, c("Melinaea_mothone",  "Hypothyris_anastasia")) # Ithomiini
  ),
  age.min = c(69.4,  6.39,  22.88),
  age.max = c(111.5, 13.77, 25.40),
  soft.bounds = TRUE
)

timetree <- chronos(tree, calibration = cal, lambda = 1,
                    model = "relaxed", quiet = TRUE)

root_age <- max(node.depth.edgelength(timetree))
message("  Root age: ", round(root_age, 1), " Ma")

# 3. Tip metadata
tip_data <- data.frame(
  label = timetree$tip.label,
  stringsAsFactors = FALSE
) %>%
  mutate(
    display_name = gsub("_", " ", label),
    colour       = SP_COLOURS[label],
    tribe = case_when(
      grepl("^Heliconius", label) ~ "Heliconiini",
      TRUE                        ~ "Ithomiini"
    )
  )

# 4. Node age labels (TimeTree high-confidence nodes only)
heli_node <- MRCA(timetree, c("Heliconius_erato",  "Heliconius_numata"))
itho_node <- MRCA(timetree, c("Melinaea_mothone",  "Hypothyris_anastasia"))
root_node <- MRCA(timetree, c("Heliconius_erato",  "Melinaea_mothone"))

root_node_label <- data.frame(
  node  = as.integer(c(root_node)),
  label = c("77\nMa"),
  stringsAsFactors = FALSE
)

inner_node_labels <- data.frame(
  node  = as.integer(c(heli_node, itho_node)),
  label = c("12\nMa", "24\nMa"),
  stringsAsFactors = FALSE
)

# 5. Build plot
message("Building plot...")

p <- ggtree(timetree, colour = "grey25", size = 0.7, layout = "rectangular") %<+% tip_data

# Extract coordinates for node label positioning
tree_data   <- p$data

root_node_coords  <- root_node_label %>%
  left_join(tree_data %>% select(node, x, y), by = "node")

inner_node_coords <- inner_node_labels %>%
  left_join(tree_data %>% select(node, x, y), by = "node")

p <- p +

  # Tree elements
  geom_tippoint(aes(colour = label), size = 3.5, show.legend = FALSE) +
  geom_tiplab(aes(label = display_name, colour = label),
              fontface = "italic", size = 6.0, offset = 1.5,
              show.legend = FALSE) +
  geom_nodepoint(colour = "grey40", size = 1.8, alpha = 0.7) +

  # Root node label as a plain ggplot2 label
  geom_label(
    data          = root_node_coords,
    aes(x = x, y = y, label = label),
    size          = 3.0,
    hjust         = -0.25,
    vjust         = 0.5,
    colour        = "black",
    fill          = "grey40",
    label.padding = unit(0.25, "lines"),
    linewidth     = 0.3,
    alpha         = 0.55,
    inherit.aes   = FALSE
  ) +

  # Inner node labels as plain ggplot2 labels
  geom_label(
    data          = inner_node_coords,
    aes(x = x, y = y, label = label),
    size          = 3.0,
    hjust         = 0.5,
    vjust         = 0.5,
    colour        = "black",
    fill          = "grey60",
    label.padding = unit(0.25, "lines"),
    linewidth     = 0.3,
    alpha         = 0.80,
    inherit.aes   = FALSE
  ) +

  # node colors
  scale_colour_manual(
    name   = "Species",
    values = SP_COLOURS,
    guide  = "none"
  ) +

  theme_tree2() +
  theme(
    legend.position = "none",
    axis.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.caption    = element_text(size = 7, colour = "grey50",
                                   hjust = 0, face = "italic"),
    plot.margin     = margin(25, 10, 10, 10)
  )

# add photo labels
photo_df <- p$data %>%
  filter(isTip) %>%
  select(label, x, y) %>%
  mutate(
    image = case_when(
      label == "Heliconius_erato"        ~ "photos/heliconius_erato.png",
      label == "Heliconius_melpomene"    ~ "photos/heliconius_melpomene.png",
      label == "Heliconius_numata"       ~ "photos/heliconius_numata.png",
      label == "Heliconius_sara"         ~ "photos/heliconius_sara.png",
      label == "Mechanitis_lysimnia"     ~ "photos/mechanitis_lysimnia.png",
      label == "Mechanitis_polymnia"     ~ "photos/mechanitis_polymnia.png",
      label == "Mechanitis_messenoides"  ~ "photos/mechanitis_messenoides.png",
      label == "Melinaea_mothone"        ~ "photos/melinaea_mothone.png",
      label == "Melinaea_menophilus"     ~ "photos/melinaea_menophilus.png",
      label == "Ithomia_salapia"         ~ "photos/ithomia_salapia.png",
      label == "Hypothyris_anastasia"    ~ "photos/hypothyris_anastasia.png",
      TRUE ~ NA_character_
    )
  )

# Add to plot
p <- p + geom_image(data = photo_df,
                    aes(x = x + 43.5, y = y, image = image),
                    size = 0.09, inherit.aes = FALSE)

# 7. Save
ggsave(
  file.path(out_dir, "species_tree_timetree.png"),
  p, width = 10, height = 8, dpi = 300, bg = "white"
)

message("Saved: species_tree_timetree.png")

write.tree(timetree, file.path(out_dir, "SpeciesTree_timetree.nwk"))
message("Saved: SpeciesTree_timetree.nwk")

message("Done. 🌳")