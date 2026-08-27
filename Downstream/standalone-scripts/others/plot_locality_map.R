#!/usr/bin/env Rscript
# plot_locality_map.R
#
# Plots sample localities on a shaded-relief topographic map
#
# Usage: Rscript plot_locality_map.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(marmap)
})

# Config
SAMPLES_CSV <- "samples_map.csv"
OUT_DIR     <- "output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read localities
samples <- read.csv(SAMPLES_CSV, stringsAsFactors = FALSE) %>%
  mutate(
    species = Species
  )

# Report any species names that failed to match the lookup table
unmatched <- samples %>% filter(is.na(species)) %>% distinct(Species)
if (nrow(unmatched) > 0) {
  message("WARNING: unrecognised species names (dropped): ",
          paste(unmatched$Species, collapse = ", "))
}

samples <- samples %>%
  filter(!is.na(species), !is.na(Latitude), !is.na(Longitude)) %>%
  mutate(
    species     = factor(species,     levels = SPECIES),
    species_lab = factor(SPECIES_LABELS[species], levels = SPECIES_LABELS[SPECIES])
  )

message("Sample counts per species:")
print(samples %>% count(species_lab))

# Basemap
buffer     <- 3
resolution <- 2

xlim <- range(samples$Longitude) + c(-buffer, buffer)
ylim <- range(samples$Latitude)  + c(-buffer, buffer)

bathy <- getNOAA.bathy(
  lon1 = xlim[1], lon2 = xlim[2],
  lat1 = ylim[1], lat2 = ylim[2],
  resolution = resolution, keep = TRUE
)

bathy_df <- fortify.bathy(bathy)  # columns: x, y, z (elevation in m)

# Only land needs to be drawn
land_df <- bathy_df %>% filter(z >= 0)

# Land colours
LAND_PAL <- c("#4c7a3d", "#8fae4e", "#b5d373", "#d9c273", "#b6813f", "#7a4a2b", "#f2f2f2")

# Shared theme
base_map_theme <- theme_map() +
  theme(
    panel.background  = element_rect(fill = "white", colour = NA),
    legend.position   = c(0.75, 0.75),
  )

# All specimens are shown as black
p1 <- ggplot() +
  geom_raster(data = land_df, aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(colours = LAND_PAL, name = "Elevation (m)", ) +
  geom_contour(
    data = bathy_df, aes(x = x, y = y, z = z),
    breaks = 0, colour = "grey20", linewidth = 0.25
  ) +
  geom_point(
    data = samples, aes(x = Longitude, y = Latitude),
    colour = "black", size = 3, alpha = 0.35
  ) +
  scale_alpha_continuous(range = c(0, 0.5), guide = "none") +
  coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = NULL, x = NULL, y = NULL) +
  coord_sf(crs = 4326) +  # WGS84
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  theme(
    axis.label = element_blank,
    axis.label.y = element_blank,
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

p1

ggsave(
  file.path(OUT_DIR, "species_distribution.pdf"),
  p1, width = 10, height = 9, bg = "white",
  limitsize = FALSE, device = cairo_pdf
)

message("Saved: ", file.path(OUT_DIR, "locality_map_combined.pdf"))

write.csv(samples, file.path(OUT_DIR, "localities_combined.csv"), row.names = FALSE)