# Figure 3E — leaf temperature depression of survivors vs bare substrate at 60 °C
# Inputs: data/fig3e_tidy_temps.csv                     (survival status per cell)
#         data/fig3e_tidy_counts_soils.csv              (substrate cells)
#         data/fig3e_radiuses_center_combined_all.csv   (IR camera temperatures)
#         data/fig3e_hobo.csv                           (chamber air T, RH, light)

library(dplyr)
library(ggplot2)
library(lmerTest)

TARGET_DATE <- "2024-10-29"     # 60 °C day
WIN_LO <- 945
WIN_HI <- 1610
TZ     <- "America/Los_Angeles"

# 0 = survivor, 1 = substrate
pal <- c("0" = "#457B9D", "1" = "#DAA520")

# --- plants and substrate cells ----------------------------------------------
plants <- read.csv("data/fig3e_tidy_temps.csv")[-1] %>%
  select(id, cell, tray, row, col, individual, fustat) %>%
  mutate(id = as.character(id), fustat = 0L)     # survivors

soils <- read.csv("data/fig3e_tidy_counts_soils.csv")[-1] %>%
  select(id, cell, tray, row, col, individual) %>%
  mutate(id = as.character(id), fustat = 1L)     # substrate

cells <- bind_rows(plants, soils)

# --- leaf and substrate temperatures -----------------------------------------
ir <- read.csv("data/fig3e_radiuses_center_combined_all.csv") %>%
  inner_join(cells, by = "cell")

t <- as.POSIXct(sub("^Record_", "", ir$csv), format = "%Y-%m-%d_%H-%M-%S", tz = TZ)
ir <- ir %>%
  mutate(csv_date     = as.Date(t, tz = TZ),
         csv_time     = format(t, "%H:%M:%S", tz = TZ),
         numeric_time = as.integer(format(t, "%H%M", tz = TZ)),
         cam_time     = as.POSIXct(paste(csv_date,
                                         format(t, "%H:%M:00", tz = TZ)), tz = TZ))

# --- chamber conditions ------------------------------------------------------
hobo <- read.csv("data/fig3e_hobo.csv", skip = 1, header = TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8-BOM")[, 1:5]
names(hobo) <- c("row", "datetime", "tempC", "rh", "light")
hobo$hobo_time <- as.POSIXct(hobo$datetime, format = "%y-%m-%d %H:%M:%S", tz = TZ)

# match each camera frame to the nearest logger reading within 3 min
idx <- sapply(ir$cam_time, function(x) {
  d <- abs(difftime(hobo$hobo_time, x, units = "mins"))
  if (all(is.na(d))) NA else which.min(d)
})

d <- cbind(ir, hobo[idx, c("tempC", "rh", "light", "hobo_time")])
d <- d[abs(as.numeric(difftime(d$hobo_time, d$cam_time, units = "mins"))) <= 3, ]

# --- radiative correction of the air sensor ----------------------------------
# The logger reads high under the lamps. The offset is the rise from pre-dawn
# (21:00-22:00) to lights-on (06:00-07:00) on the 55 °C day.
window_mean <- function(date, lo, hi) {
  x <- subset(hobo, as.Date(hobo_time) == as.Date(date))
  h <- as.numeric(format(x$hobo_time, "%H")) +
    as.numeric(format(x$hobo_time, "%M")) / 60
  mean(x$tempC[h >= lo & h < hi], na.rm = TRUE)
}

avg_dev <- window_mean("2024-10-21", 6, 7) - window_mean("2024-10-21", 21, 22)

d$tempC_adj <- d$tempC - avg_dev
d$deltaT    <- d$area_mean_3 - d$tempC_adj

# --- analysis window ---------------------------------------------------------
w <- d %>%
  filter(csv_date == as.Date(TARGET_DATE),
         numeric_time >= WIN_LO, numeric_time <= WIN_HI) %>%
  mutate(time_posix = as.POSIXct(substr(csv_time, 1, 5), format = "%H:%M", tz = "UTC"),
         grp        = factor(fustat, levels = c(1, 0)))

offset <- mean(w$area_mean_3 - w$deltaT, na.rm = TRUE)   # ΔT -> Tleaf axis

# --- left panel: time series -------------------------------------------------
model <- lmer(deltaT ~ grp * scale(as.numeric(time_posix)) + (1 | cell),
              data = w)
print(anova(model))

p_group <- anova(model)["grp", "Pr(>F)"]
cap <- paste0("p(group) = ", formatC(p_group, format = "e", digits = 2))

agg <- w %>%
  group_by(grp, time_posix) %>%
  summarise(mean = mean(deltaT, na.rm = TRUE),
            sd   = sd(deltaT, na.rm = TRUE),
            n    = sum(!is.na(deltaT)), .groups = "drop") %>%
  mutate(ci    = 1.96 * sd / sqrt(pmax(n, 1)),
         lower = mean - ci, upper = mean + ci)

p_series <- ggplot(agg, aes(time_posix, mean, color = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
  scale_y_continuous(
    name     = expression(Delta * italic(T) ~ "(°C)"),
    sec.axis = sec_axis(~ . + offset, name = expression(italic(T)[leaf] ~ "(°C)"))
  ) +
  scale_color_manual(values = pal, labels = c("1" = "Substrate", "0" = "Survivors")) +
  scale_fill_manual(values  = pal, labels = c("1" = "Substrate", "0" = "Survivors")) +
  annotate("text", x = max(agg$time_posix), y = -Inf, label = cap,
           hjust = 1.05, vjust = -1, size = 4, family = "Helvetica",
           color = "grey30") +
  labs(x = "Time of Day", color = "Group", fill = "Group") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 13),
    legend.text        = element_text(size = 12),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    axis.text.x        = element_text(size = 8),
    axis.text.y        = element_text(size = 12),
    axis.title         = element_text(size = 14),
    axis.title.y.left  = element_text(margin = margin(r = -5)),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    panel.background   = element_rect(fill = "transparent", color = NA),
    plot.background    = element_rect(fill = "transparent", color = NA)
  )

pdf("fig3e_left_deltaT_series.pdf", width = 5, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_series)
dev.off()

# --- right panel: per-cell extremes ------------------------------------------
s <- w %>%
  group_by(cell, grp) %>%
  summarise(min = min(deltaT, na.rm = TRUE),
            max = max(deltaT, na.rm = TRUE), .groups = "drop")

long <- bind_rows(
  data.frame(grp = s$grp, value = s$min, stat = "min"),
  data.frame(grp = s$grp, value = s$max, stat = "max")
)

p_box <- ggplot(long, aes(grp, value, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_fill_manual(values = pal, labels = c("1" = "Substrate", "0" = "Survivors")) +
  scale_y_continuous(
    breaks   = seq(-10, 10, by = 3),
    limits   = c(-11, 6),
    labels   = function(x) gsub("-", "\u2212", x),
    name     = expression(Delta * italic(T) ~ "(°C)"),
    sec.axis = sec_axis(~ . + offset, name = expression(italic(T)[leaf] ~ "(°C)"),
                        breaks = seq(round(-10 + offset), round(10 + offset), by = 3))
  ) +
  facet_wrap(~factor(stat, levels = c("max", "min")), strip.position = "bottom",
             labeller = as_labeller(c("min" = "Abs Max", "max" = "Abs Min"))) +
  labs(x = NULL, fill = "Group") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position    = "none",
    strip.background   = element_blank(),
    strip.placement    = "outside",
    strip.text         = element_text(size = 12),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 12),
    axis.title         = element_text(size = 14),
    axis.title.y.left  = element_text(margin = margin(r = -5)),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    panel.background   = element_rect(fill = "transparent", color = NA),
    plot.background    = element_rect(fill = "transparent", color = NA)
  )

pdf("fig3e_right_deltaT_extremes.pdf", width = 2.65, height = 2.7,
    useDingbats = FALSE, bg = "transparent")
print(p_box)
dev.off()