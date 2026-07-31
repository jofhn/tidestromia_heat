# Figure 1A — T. oblongifolia occurrences over summer maximum temperature
# Occurrence records and climate rasters are fetched at runtime; no input files.

library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(geodata)
library(rgbif)
library(rinat)

# --- occurrences -------------------------------------------------------------
# Mahalanobis outlier removal, applied separately to each source
drop_outliers <- function(df, p) {
  m <- as.matrix(df[, c("decimalLongitude", "decimalLatitude")])
  d <- mahalanobis(m, colMeans(m), cov(m))
  df[d <= qchisq(p, df = 2), ]
}

gbif <- occ_search(scientificName = "Tidestromia oblongifolia",
                   hasCoordinate = TRUE, limit = 10000)$data %>%
  select(decimalLongitude, decimalLatitude) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  drop_outliers(0.95) %>%
  mutate(source = "GBIF")

inat <- get_inat_obs(taxon_name = "Tidestromia oblongifolia",
                     quality = "research", maxresults = 10000) %>%
  select(decimalLongitude = longitude, decimalLatitude = latitude) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  drop_outliers(0.999) %>%
  mutate(source = "iNaturalist")

occ_sf <- st_as_sf(rbind(gbif, inat),
                   coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# --- extent, basemap, climate ------------------------------------------------
bb <- st_bbox(st_as_sf(gbif, coords = c("decimalLongitude", "decimalLatitude"),
                       crs = 4326))
xmin <- bb["xmin"] - 1; xmax <- bb["xmax"] + 1
ymin <- bb["ymin"] - 1; ymax <- bb["ymax"] + 1

usa    <- ne_states(country = "united states of america", returnclass = "sf")
mexico <- ne_states(country = "mexico", returnclass = "sf")

tmax <- worldclim_global(var = "tmax", res = 2.5, path = ".")
summer_crop <- crop(mean(tmax[[c(6, 7, 8)]]), ext(xmin, xmax, ymin, ymax))

landmarks <- data.frame(
  name = c("Death Valley", "Las Vegas", "Anza-Borrego"),
  lon  = c(-116.87, -115.14, -116.24),
  lat  = c(36.46, 36.17, 33.09)
)

# --- plot --------------------------------------------------------------------
p <- ggplot() +
  geom_spatraster(data = summer_crop) +
  scale_fill_gradientn(
    colors = c("#2255cc", "#2299dd", "#33bbaa", "#99bb22",
               "#ddcc22", "#dd7722", "#cc0000"),
    limits = c(10, 45),
    name = expression(italic(T)[max] ~ "(°C)"),
    na.value = "transparent"
  ) +
  geom_sf(data = usa,    fill = NA, color = "lightgrey", linewidth = 0.3) +
  geom_sf(data = mexico, fill = NA, color = "lightgrey", linewidth = 0.3) +
  geom_sf(data = occ_sf, aes(color = source), size = 1, alpha = 0.5) +
  scale_color_manual(values = c(GBIF = "black", iNaturalist = "darkgrey"),
                     name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  geom_point(data = landmarks, aes(lon, lat), color = "white", size = 2) +
  geom_text(data = landmarks, aes(lon, lat, label = name),
            color = "white", size = 3.5, family = "Helvetica",
            nudge_x = c(0.2, 1.8, 1.4), nudge_y = c(0.55, 0.3, -0.5)) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           expand = FALSE, clip = "on") +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.4, "cm"),
    axis.title       = element_text(size = 12),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 8, angle = 30, hjust = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("fig1a_range_map.pdf", width = 5, height = 3.12, useDingbats = FALSE,
    bg = "transparent")
print(p)
dev.off()