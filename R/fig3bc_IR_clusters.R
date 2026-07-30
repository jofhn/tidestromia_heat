# Figure 3B-C — leaf temperature depression by cluster at 60 °C
# Inputs: data/fig3abc_tidy_temps.csv          (cell metadata, clusters)
#         data/fig3bcd_radiuses_leaves_combined_all.csv    (IR camera leaf temperatures)
#         data/fig3bc_hobo.csv (chamber air T, RH, light)

library(dplyr)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(ggtext)
library(glue)

TARGET_DATE <- "2025-04-16"      # 60 °C day
WIN_LO      <- 1035              # analysis window, HHMM
WIN_HI      <- 1500
TZ          <- "America/Los_Angeles"

pal <- c("1" = "#E63946", "2" = "#F4A261", "3" = "#457B9D")

# --- leaf temperatures -------------------------------------------------------
meta <- read.csv("data/fig3abc_tidy_temps.csv")[-1] %>%
  select(id, cluster, tray, row, col, cell, si, site, individual,
         family, where, time_to_four)

comprehensive <- read.csv("data/fig3bcd_radiuses_leaves_combined_all.csv") %>%
  merge(meta, by = "cell") %>%
  mutate(cluster = as.factor(cluster))

t <- as.POSIXct(sub("^Record_", "", comprehensive$csv),
                format = "%Y-%m-%d_%H-%M-%S", tz = TZ)
comprehensive <- comprehensive %>%
  mutate(csv_date     = as.Date(t, tz = TZ),
         csv_time     = format(t, "%H:%M:%S", tz = TZ),
         numeric_time = as.integer(format(t, "%H%M", tz = TZ)),
         cam_time     = as.POSIXct(paste(csv_date,
                                         format(t, "%H:%M:00", tz = TZ)), tz = TZ))

# --- chamber conditions ------------------------------------------------------
hobo <- read.csv("data/fig3bc_hobo.csv", skip = 1,
                 header = TRUE, check.names = FALSE,
                 fileEncoding = "UTF-8-BOM")[, 1:5]
names(hobo) <- c("row", "datetime", "tempC", "rh", "light")
hobo$hobo_time <- as.POSIXct(hobo$datetime, format = "%y-%m-%d %H:%M:%S", tz = TZ)

# match each camera frame to the nearest logger reading within 3 min
idx <- sapply(comprehensive$cam_time, function(x) {
  d <- abs(difftime(hobo$hobo_time, x, units = "mins"))
  if (all(is.na(d))) NA else which.min(d)
})

comprehensive <- cbind(comprehensive, hobo[idx, c("tempC", "rh", "light", "hobo_time")])
comprehensive <- comprehensive[
  abs(as.numeric(difftime(comprehensive$hobo_time,
                          comprehensive$cam_time, units = "mins"))) <= 3, ]

# --- radiative correction of the air sensor ----------------------------------
# The logger reads high when the lamps are on. The offset is estimated as the
# rise from pre-dawn (21:00-22:00 previous night) to lights-on (06:00-07:00),
# averaged over the two hottest days, and applied only while lights are on.
window_mean <- function(date, lo, hi) {
  d <- subset(hobo, as.Date(hobo_time) == as.Date(date))
  h <- as.numeric(format(d$hobo_time, "%H")) +
    as.numeric(format(d$hobo_time, "%M")) / 60
  mean(d$tempC[h >= lo & h < hi], na.rm = TRUE)
}

avg_dev <- mean(c(
  window_mean("2025-04-15", 6, 7) - window_mean("2025-04-15", 21, 22),
  window_mean("2025-04-16", 6, 7) - window_mean("2025-04-16", 21, 22)
))

lights_on <- comprehensive$light > 1000      # dark floor ~11.8 lux, lit ~32280
comprehensive$tempC_adj <- comprehensive$tempC
comprehensive$tempC_adj[lights_on] <- comprehensive$tempC[lights_on] - avg_dev

comprehensive$deltaT <- comprehensive$area_mean_3 - comprehensive$tempC_adj

# --- model -------------------------------------------------------------------
dd <- comprehensive %>%
  filter(csv_date == as.Date(TARGET_DATE),
         numeric_time >= WIN_LO, numeric_time <= WIN_HI,
         cluster != "soil") %>%
  mutate(time_posix  = as.POSIXct(substr(csv_time, 1, 5), format = "%H:%M", tz = "UTC"),
         time_scaled = scale(as.numeric(time_posix)),
         grp_fac     = droplevels(as.factor(cluster)),
         cell_fac    = as.factor(cell))

model <- lmer(deltaT ~ grp_fac * time_scaled + (1 | cell_fac), data = dd)
print(anova(model))

p_cluster <- anova(model)["grp_fac", "Pr(>F)"]
cap <- paste0("p(cluster) = ", formatC(p_cluster, format = "e", digits = 2))

# --- 3B: time series ---------------------------------------------------------
agg <- dd %>%
  group_by(grp = as.character(cluster), time_posix) %>%
  summarise(mean = mean(deltaT, na.rm = TRUE),
            sd   = sd(deltaT, na.rm = TRUE),
            n    = sum(!is.na(deltaT)), .groups = "drop") %>%
  mutate(ci    = 1.96 * sd / sqrt(pmax(n, 1)),
         lower = mean - ci, upper = mean + ci)

offset <- mean(dd$area_mean_3 - dd$deltaT, na.rm = TRUE)   # ΔT -> Tleaf axis

p_series <- ggplot(agg, aes(time_posix, mean, color = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hours",
                   date_minor_breaks = "30 mins") +
  scale_y_continuous(
    name     = expression(Delta * italic(T) ~ "(°C)"),
    sec.axis = sec_axis(~ . + offset, name = expression(italic(T)[leaf] ~ "(°C)"))
  ) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  annotate("text", x = max(agg$time_posix), y = -Inf, label = cap,
           hjust = 1.05, vjust = -1, size = 4, family = "Helvetica",
           color = "grey30") +
  labs(x = "Time of Day") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position    = "none",
    axis.text          = element_text(size = 12),
    axis.title         = element_text(size = 14),
    axis.title.y.left  = element_text(margin = margin(r = -5)),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    panel.background   = element_rect(fill = "transparent", color = NA),
    plot.background    = element_rect(fill = "transparent", color = NA)
  )

pdf("fig3b_deltaT_series.pdf", width = 4, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_series)
dev.off()

# --- 3C: pairwise cluster contrasts ------------------------------------------
pw <- as.data.frame(contrast(emmeans(model, ~ grp_fac), method = list(
  "3 - 2" = c( 0, -1,  1),
  "3 - 1" = c(-1,  0,  1),
  "2 - 1" = c(-1,  1,  0)
), adjust = "none")) %>%
  mutate(
    conf.low  = estimate - 1.96 * SE,
    conf.high = estimate + 1.96 * SE,
    sig = case_when(p.value < 0.001 ~ "***",
                    p.value < 0.01  ~ "**",
                    p.value < 0.05  ~ "*",
                    TRUE            ~ "ns"),
    contrast_label = case_when(
      contrast == "3 - 2" ~ glue("<span style='color:{pal['3']}'>**3**</span> - <span style='color:{pal['2']}'>**2**</span>"),
      contrast == "3 - 1" ~ glue("<span style='color:{pal['3']}'>**3**</span> - <span style='color:{pal['1']}'>**1**</span>"),
      contrast == "2 - 1" ~ glue("<span style='color:{pal['2']}'>**2**</span> - <span style='color:{pal['1']}'>**1**</span>")
    )
  )

p_coef <- ggplot(pw, aes(estimate, contrast_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.15) +
  geom_point(size = 3, shape = 16) +
  geom_text(aes(x = conf.high, label = paste0(" ", sig)),
            hjust = 0, size = 4.5, family = "Helvetica") +
  coord_cartesian(clip = "off") +
  labs(x = expression("Mean difference in " * Delta * italic(T) ~ "(°C)"),
       y = "Cluster") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position   = "none",
    axis.text.x       = element_text(size = 12),
    axis.text.y       = element_markdown(family = "Helvetica", size = 12),
    axis.title        = element_text(size = 14),
    axis.title.y.left = element_text(margin = margin(r = 5)),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    plot.margin       = margin(l = 20)
  )

pdf("fig3c_cluster_contrasts.pdf", width = 3.525, height = 2.65,
    useDingbats = FALSE, bg = "transparent")
print(p_coef)
dev.off()
