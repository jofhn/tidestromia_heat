# Figure 1B-D — collection sites and ERA5 summer climate by region
# Inputs: data/collection_regions.csv
#         data/era5_land_temperature_2m.csv  (436 MB; in Zenodo deposit or regenerate with R/fetch_era_land5.R)

library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(geodata)
library(rnaturalearth)
library(ggnewscale)
library(lmerTest)
library(emmeans)
library(multcompView)

pal <- c(DV = "#EE6677", perimeter = "#4477AA", outside = "#228833")
regions <- c("DV", "perimeter", "outside")

# --- data --------------------------------------------------------------------
csv <- "data/climates/era5_land_temperature_2m.csv"
if (!file.exists(csv)) {
  stop("Missing ", csv, " — download it from the Zenodo deposit, or regenerate ",
       "it with R/fetch_era_land5.R")
}

weather <- read.csv(csv) %>%
  select(-any_of("X")) %>%
  mutate(
    time  = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    month = factor(format(time, "%m"), levels = c("06", "07", "08"),
                   labels = c("June", "July", "August")),
    hour  = as.numeric(format(time, "%H")),
    date  = as.Date(time)
  ) %>%
  filter(where %in% regions)

sites_sf <- st_as_sf(weather, coords = c("long", "lat"), crs = 4326)

# --- 1B: terrain map ---------------------------------------------------------
xmin <- -118.5; xmax <- -111.5; ymin <- 31.5; ymax <- 38.5   # see note below

elev <- terra::merge(elevation_30s(country = "USA", path = "geodata"),
                     elevation_30s(country = "MEX", path = "geodata")) %>%
  crop(ext(xmin - 0.5, xmax + 0.5, ymin - 0.5, ymax + 0.5))

hill <- shade(terrain(elev, "slope",  unit = "radians"),
              terrain(elev, "aspect", unit = "radians"),
              angle = 40, direction = 315)

p_terrain <- ggplot() +
  geom_spatraster(data = hill) +
  scale_fill_distiller(palette = "Greys", direction = -1,
                       na.value = "transparent", guide = "none") +
  geom_sf(data = ne_states(country = "united states of america",
                           returnclass = "sf"),
          fill = NA, color = "lightgrey", linewidth = 0.3) +
  geom_sf(data = ne_countries(country = "mexico", scale = "medium",
                              returnclass = "sf"),
          fill = NA, color = "lightgrey", linewidth = 0.3) +
  new_scale_color() +
  geom_sf(data = sites_sf, aes(color = where), shape = 16, size = 3, alpha = 0.9) +
  scale_color_manual(values = pal[regions], breaks = regions,
                     labels = c("Death Valley", "Perimeter", "Outside"),
                     name = "Collection Region") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 5))) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           expand = FALSE, crs = 4326) +
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
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# hillshade is a raster, so export via PNG and convert
ragg::agg_png(tmp <- tempfile(fileext = ".png"), width = 4.25, height = 3.12,
              units = "in", res = 300, background = "transparent")
print(p_terrain); dev.off()
magick::image_write(magick::image_read(tmp), "fig1b_terrain.pdf",
                    format = "pdf", density = 300)

# --- 1C: mean hourly temperature, daytime hours ------------------------------
model_day <- lmer(
  temperature_2m ~ where * month + hour + (1 | site),
  data = weather %>%
    filter(hour >= 11, hour <= 21) %>%
    mutate(where = factor(where, levels = regions),
           month = factor(month, levels = c("July", "June", "August")))
)
print(anova(model_day))
print(pairs(emmeans(model_day, ~ where | month), adjust = "tukey", infer = TRUE))

temp_day <- weather %>%
  filter(hour >= 11, hour <= 21) %>%
  group_by(where, month, hour) %>%
  summarise(avg_temp = mean(temperature_2m, na.rm = TRUE),
            sd_temp  = sd(temperature_2m, na.rm = TRUE), .groups = "drop")

p_curves <- ggplot(temp_day, aes(hour, avg_temp, color = where, group = where)) +
  geom_ribbon(aes(ymin = avg_temp - sd_temp, ymax = avg_temp + sd_temp,
                  fill = where), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~month) +
  scale_x_continuous(breaks = seq(0, 23, 2)) +
  scale_color_manual(values = pal, limits = regions) +
  scale_fill_manual(values = pal, limits = regions) +
  scale_y_continuous(limits = c(24, 47), breaks = seq(25, 45, 5)) +
  labs(x = "Time of Day", y = expression(italic(T)[avg] ~ "(°C)"),
       color = NULL, fill = NULL) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 13),
    legend.position  = "none",
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 12),
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 16),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("fig1c_temp_curves.pdf", width = 4.75, height = 2.8,
    useDingbats = FALSE, bg = "transparent")
print(p_curves); dev.off()

# --- 1D: daily maximum temperature -------------------------------------------
daily_tmax <- weather %>%
  group_by(site, where, month, date) %>%
  summarise(daily_tmax = max(temperature_2m, na.rm = TRUE), .groups = "drop")

model_tmax <- lmer(daily_tmax ~ where * month + (1 | site), data = daily_tmax)
print(anova(model_tmax))

emm_df <- as.data.frame(pairs(emmeans(model_tmax, ~ where | month),
                              adjust = "tukey"))

letters_tmax <- lapply(levels(daily_tmax$month), function(m) {
  sub <- emm_df[emm_df$month == m, ]
  p <- setNames(sub$p.value, gsub(" - ", "-", trimws(as.character(sub$contrast))))
  l <- multcompLetters(p)$Letters
  data.frame(where = names(l), Letters = l, month = m)
}) %>% bind_rows()

label_positions <- daily_tmax %>%
  group_by(where, month) %>%
  summarise(y_max = quantile(daily_tmax, 0.75) + 1.5 * IQR(daily_tmax),
            .groups = "drop") %>%
  left_join(letters_tmax, by = c("where", "month"))

p_tmax <- ggplot(daily_tmax, aes(factor(where, levels = regions),
                                 daily_tmax, fill = where)) +
  geom_boxplot(outlier.shape = NA, width = 0.5,
               color = "grey30", linewidth = 0.5) +
  facet_wrap(~month) +
  scale_fill_manual(values = pal, limits = regions) +
  scale_x_discrete(limits = regions,
                   labels = c("DV", "Perimeter", "Outside")) +
  scale_y_continuous(limits = c(20, 56), breaks = seq(20, 60, 5)) +
  coord_cartesian(clip = "off") +
  geom_text(data = label_positions, aes(x = where, label = Letters),
            y = Inf, vjust = 0.5, size = 5,
            fontface = "bold", family = "Helvetica") +
  labs(x = "Population", y = expression(italic(T)[max] ~ "(°C)")) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y      = element_text(size = 12),
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 16),
    strip.background = element_blank(),
    strip.text       = element_text(size = 14, vjust = 3),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("fig1d_tmax_box.pdf", width = 3.5, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_tmax); dev.off()