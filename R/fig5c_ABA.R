# Figure 5C — leaf temperature depression under ABA at 60 °C
# Inputs: data/fig5c_tidy_counts.csv                  (cell metadata)
#         data/fig5c_treatments.csv                   (ABA / mock / untreated)
#         data/fig5c_radiuses_leaves_combined_all.csv (IR camera temperatures)
#         data/fig5c_hobo.csv                         (chamber air T, RH, light)

library(dplyr)
library(ggplot2)
library(lmerTest)
library(multcompView)

TARGET_DATE <- "2025-10-06"   # 60 °C day
WIN_LO <- 1115                # one hour spanning peak leaf cooling
WIN_HI <- 1215
TZ     <- "America/Los_Angeles"

pal <- c(aba = "#E63946", mock = "#8B7BAD", untreated = "#A0785A")
lab <- c(aba = "ABA", mock = "Mock", untreated = "Untreated")

# --- plants and treatments ---------------------------------------------------
meta <- read.csv("data/fig5c_tidy_counts.csv")[-1]
treatments <- read.csv("data/fig5c_treatments.csv")

cells <- merge(meta, treatments, by = c("tray", "row", "col", "cell", "id")) %>%
  rename(treatment = round) %>%
  select(id, cell, tray, row, col, individual, treatment) %>%
  filter(treatment %in% names(pal))     # drops the "fca" arm and substrate cells

# --- leaf temperatures -------------------------------------------------------
ir <- read.csv("data/fig5c_radiuses_leaves_combined_all.csv") %>%
  inner_join(cells, by = "cell")

t <- as.POSIXct(sub("^Record_", "", ir$csv), format = "%Y-%m-%d_%H-%M-%S", tz = TZ)
ir <- ir %>%
  mutate(csv_date     = as.Date(t, tz = TZ),
         csv_time     = format(t, "%H:%M:%S", tz = TZ),
         numeric_time = as.integer(format(t, "%H%M", tz = TZ)),
         cam_time     = as.POSIXct(paste(csv_date,
                                         format(t, "%H:%M:00", tz = TZ)), tz = TZ))

# --- chamber conditions ------------------------------------------------------
hobo <- read.csv("data/fig5c_hobo.csv", skip = 1, header = TRUE,
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
# (21:30-22:30) to lights-on (06:30-07:30) on the 55 °C day.
window_mean <- function(date, lo, hi) {
  x <- subset(hobo, as.Date(hobo_time) == as.Date(date))
  h <- as.numeric(format(x$hobo_time, "%H")) +
    as.numeric(format(x$hobo_time, "%M")) / 60
  mean(x$tempC[h >= lo & h < hi], na.rm = TRUE)
}

avg_dev <- window_mean("2025-10-05", 6.5, 7.5) - window_mean("2025-10-05", 21.5, 22.5)

d$tempC_adj <- d$tempC - avg_dev
d$deltaT    <- d$area_mean_3 - d$tempC_adj

# --- analysis window ---------------------------------------------------------
w <- d %>%
  filter(csv_date == as.Date(TARGET_DATE),
         numeric_time >= WIN_LO, numeric_time <= WIN_HI) %>%
  mutate(time_posix = as.POSIXct(substr(csv_time, 1, 5), format = "%H:%M", tz = "UTC"),
         treatment  = factor(treatment, levels = names(pal)))

offset <- mean(w$area_mean_3 - w$deltaT, na.rm = TRUE)   # ΔT -> Tleaf axis

# --- left panel: time series -------------------------------------------------
model <- lmer(deltaT ~ treatment * scale(as.numeric(time_posix)) + (1 | cell),
              data = w)
print(anova(model))

p_treat <- anova(model)["treatment", "Pr(>F)"]
cap <- paste0("p(treatment) = ",
              if (p_treat < 2.2e-16) "< 2.2e-16"
              else formatC(p_treat, format = "e", digits = 2))

agg <- w %>%
  group_by(treatment, time_posix) %>%
  summarise(mean = mean(deltaT, na.rm = TRUE),
            sd   = sd(deltaT, na.rm = TRUE),
            n    = sum(!is.na(deltaT)), .groups = "drop") %>%
  mutate(ci    = 1.96 * sd / sqrt(pmax(n, 1)),
         lower = mean - ci, upper = mean + ci)

p_series <- ggplot(agg, aes(time_posix, mean, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "15 mins") +
  scale_y_continuous(
    breaks = seq(-12, -8, 1), limits = c(-12, -8),
    labels = function(x) gsub("-", "\u2212", x),
    name   = expression(Delta * italic(T) ~ "(°C)"),
    sec.axis = sec_axis(~ . + offset, name = expression(italic(T)[leaf] ~ "(°C)"),
                        labels = function(x) gsub("-", "\u2212", x))
  ) +
  scale_color_manual(values = pal, labels = lab) +
  scale_fill_manual(values = pal, labels = lab) +
  annotate("text", x = max(agg$time_posix), y = -Inf, label = cap,
           hjust = 0.99, vjust = -0.5, size = 3.5, family = "Helvetica",
           color = "grey30") +
  labs(x = "Time of Day", color = NULL, fill = NULL) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 12),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    legend.margin      = margin(b = -2),
    legend.box.margin  = margin(b = -2),
    axis.text.x        = element_text(size = 10),
    axis.text.y        = element_text(size = 12),
    axis.title         = element_text(size = 14),
    axis.title.y.left  = element_text(margin = margin(r = -2)),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    panel.background   = element_rect(fill = "transparent", color = NA),
    plot.background    = element_rect(fill = "transparent", color = NA)
  )

pdf("fig5c_left_aba_series.pdf", width = 3.4, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_series); dev.off()

# --- right panel: per-plant extremes -----------------------------------------
s <- w %>%
  group_by(cell, treatment) %>%
  summarise(min = min(deltaT, na.rm = TRUE),
            max = max(deltaT, na.rm = TRUE), .groups = "drop")

letters_df <- bind_rows(lapply(c("min", "max"), function(stat) {
  fit <- aov(as.formula(paste(stat, "~ treatment")), data = s)
  pv  <- TukeyHSD(fit)[[1]][, "p adj"]
  names(pv) <- rownames(TukeyHSD(fit)[[1]])
  L <- multcompLetters(pv)$Letters
  ypos <- tapply(s[[stat]], s$treatment, quantile, 0.95, na.rm = TRUE)
  data.frame(treatment = names(L), Letters = unname(L), stat = stat,
             y = ypos[names(L)] + if (stat == "min") 1 else 0.75)
}))

print(summary(aov(min ~ treatment, data = s)))
print(TukeyHSD(aov(min ~ treatment, data = s)))

long <- bind_rows(
  data.frame(treatment = s$treatment, value = s$min, stat = "min"),
  data.frame(treatment = s$treatment, value = s$max, stat = "max")
)

p_box <- ggplot(long, aes(treatment, value, fill = treatment)) +
  geom_boxplot(data = subset(long, stat == "max"), outlier.shape = NA,
               width = 0.5, color = "grey30", alpha = 0.6) +
  geom_boxplot(data = subset(long, stat == "min"), outlier.shape = NA,
               width = 0.5, color = "grey30") +
  geom_text(data = letters_df, aes(treatment, y, label = Letters),
            inherit.aes = FALSE, size = 5, fontface = "bold",
            family = "Helvetica") +
  scale_fill_manual(values = pal, labels = lab) +
  scale_y_continuous(
    breaks = seq(-13, -6, 1), limits = c(-13, -6),
    labels = function(x) gsub("-", "\u2212", x),
    name   = expression(Delta * italic(T) ~ "(°C)"),
    sec.axis = sec_axis(~ . + offset, name = expression(italic(T)[leaf] ~ "(°C)"),
                        breaks = seq(round(-13 + offset), round(-6 + offset), 1))
  ) +
  facet_wrap(~factor(stat, levels = c("max", "min")), strip.position = "bottom",
             labeller = as_labeller(c("max" = "Abs Min", "min" = "Abs Max"))) +
  labs(x = NULL, fill = "Treatment") +
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
    axis.title.y.left  = element_text(margin = margin(r = 2)),
    axis.title.y.right = element_text(margin = margin(l = 5)),
    panel.background   = element_rect(fill = "transparent", color = NA),
    plot.background    = element_rect(fill = "transparent", color = NA)
  )

pdf("fig5c_right_aba_extremes.pdf", width = 2.75, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_box); dev.off()