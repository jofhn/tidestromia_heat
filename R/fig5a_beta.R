# Figure 5A — leaf temperature depression across air temperatures, by cluster
# Inputs: data/fig3abc_tidy_temps.csv
#         data/fig3bcd_radiuses_leaves_combined_all.csv
#         data/fig3bc_hobo.csv

library(dplyr)
library(ggplot2)
library(multcompView)
library(ggh4x)

TZ <- "America/Los_Angeles"

# 51, 55 and 60 °C days, with the analysis window for each
DAYS <- list("2025-04-14" = c(1025, 1500),
             "2025-04-15" = c(1025, 1500),
             "2025-04-16" = c(1030, 1500))
DAY_LABELS <- c("1" = "51°C", "2" = "55°C", "3" = "60°C")

pal <- c("1" = "#E63946", "2" = "#F4A261", "3" = "#457B9D")

# --- data (as in fig3bc_IR_clusters.R) ---------------------------------------
meta <- read.csv("data/fig3abc_tidy_temps.csv")[-1] %>%
  select(id, cluster, cell)

ir <- read.csv("data/fig3bcd_radiuses_leaves_combined_all.csv") %>%
  merge(meta, by = "cell")

t <- as.POSIXct(sub("^Record_", "", ir$csv), format = "%Y-%m-%d_%H-%M-%S", tz = TZ)
ir <- ir %>%
  mutate(csv_date     = as.Date(t, tz = TZ),
         numeric_time = as.integer(format(t, "%H%M", tz = TZ)),
         cam_time     = as.POSIXct(paste(csv_date,
                                         format(t, "%H:%M:00", tz = TZ)), tz = TZ))

hobo <- read.csv("data/fig3bc_hobo.csv", skip = 1, header = TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8-BOM")[, 1:5]
names(hobo) <- c("row", "datetime", "tempC", "rh", "light")
hobo$hobo_time <- as.POSIXct(hobo$datetime, format = "%y-%m-%d %H:%M:%S", tz = TZ)

idx <- sapply(ir$cam_time, function(x) {
  d <- abs(difftime(hobo$hobo_time, x, units = "mins"))
  if (all(is.na(d))) NA else which.min(d)
})

d <- cbind(ir, hobo[idx, c("tempC", "light", "hobo_time")])
d <- d[abs(as.numeric(difftime(d$hobo_time, d$cam_time, units = "mins"))) <= 3, ]

# the air sensor reads high under the lamps; correct only while lights are on
window_mean <- function(date, lo, hi) {
  x <- subset(hobo, as.Date(hobo_time) == as.Date(date))
  h <- as.numeric(format(x$hobo_time, "%H")) +
    as.numeric(format(x$hobo_time, "%M")) / 60
  mean(x$tempC[h >= lo & h < hi], na.rm = TRUE)
}

avg_dev <- mean(c(
  window_mean("2025-04-15", 6, 7) - window_mean("2025-04-15", 21, 22),
  window_mean("2025-04-16", 6, 7) - window_mean("2025-04-16", 21, 22)
))

lights_on <- d$light > 1000        # dark floor ~11.8 lux, lit ~32280
d$tempC_adj <- d$tempC
d$tempC_adj[lights_on] <- d$tempC[lights_on] - avg_dev
d$deltaT <- d$area_mean_3 - d$tempC_adj

# --- per-plant extremes and mean air temperature, per day --------------------
daily <- bind_rows(lapply(seq_along(DAYS), function(i) {
  w <- DAYS[[i]]
  d %>%
    filter(csv_date == as.Date(names(DAYS)[i]),
           numeric_time >= w[1], numeric_time <= w[2],
           !is.na(deltaT), cluster != "soil") %>%
    group_by(id, cluster) %>%
    summarise(deltaT_min = min(deltaT), deltaT_max = max(deltaT),
              Tair = mean(tempC_adj), .groups = "drop") %>%
    mutate(day = i)
})) %>%
  mutate(cluster = factor(cluster, levels = names(pal)))

# --- left: extremes by air temperature ---------------------------------------
long <- bind_rows(
  daily %>% transmute(cluster, day, value = deltaT_min, stat = "min"),
  daily %>% transmute(cluster, day, value = deltaT_max, stat = "max")
) %>%
  mutate(stat = factor(stat, levels = c("max", "min")))

# Tukey letters within each day and statistic
letters_df <- bind_rows(lapply(c("min", "max"), function(st) {
  bind_rows(lapply(unique(long$day), function(dy) {
    sub <- long[long$stat == st & long$day == dy, ]
    pv  <- TukeyHSD(aov(value ~ cluster, data = sub))[[1]][, "p adj"]
    names(pv) <- rownames(TukeyHSD(aov(value ~ cluster, data = sub))[[1]])
    L <- multcompLetters(pv)$Letters
    data.frame(cluster = names(L), letter = unname(L), stat = st, day = dy)
  }))
})) %>%
  mutate(cluster = factor(cluster, levels = names(pal)),
         stat    = factor(stat, levels = c("max", "min")),
         y       = Inf)

p_minmax <- ggplot(long, aes(cluster, value, fill = cluster)) +
  geom_boxplot(data = subset(long, stat == "max"), outlier.shape = NA,
               width = 0.5, color = "grey30", alpha = 0.6) +
  geom_boxplot(data = subset(long, stat == "min"), outlier.shape = NA,
               width = 0.5, color = "grey30") +
  geom_text(data = letters_df,
            aes(cluster, y, label = letter, color = cluster),
            inherit.aes = FALSE, size = 3, fontface = "bold",
            family = "Helvetica", hjust = 0.5, vjust = 1.5) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  scale_x_discrete(labels = NULL) +
  facet_grid(stat ~ day, scales = "free_y", switch = "x",
             labeller = labeller(day = DAY_LABELS,
                                 stat = c("max" = "Abs Min", "min" = "Abs Max"))) +
  facetted_pos_scales(y = list(
    stat == "max" ~ scale_y_continuous(
      name = expression(Delta * italic(T) ~ "(°C)"), limits = c(-7, 0),
      breaks = seq(-7, 0, 2), labels = function(x) gsub("-", "\u2212", x)),
    stat == "min" ~ scale_y_continuous(
      name = expression(Delta * italic(T) ~ "(°C)"), limits = c(-13, -6),
      breaks = seq(-13, -6, 2), labels = function(x) gsub("-", "\u2212", x))
  )) +
  guides(color = "none") +
  labs(x = expression(italic(T)[air] ~ "(°C)"), fill = "Cluster") +
  theme_classic(base_family = "Helvetica", base_size = 12) +
  theme(
    strip.background  = element_blank(),
    strip.placement   = "outside",
    strip.text        = element_text(size = 14),
    legend.position   = "top",
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 12),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.key.size   = unit(0.4, "cm"),
    legend.margin     = margin(t = -5, b = -5),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 12),
    axis.title        = element_text(size = 14),
    axis.title.y.left = element_text(margin = margin(r = -2)),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA)
  )

pdf("fig5a_left_deltaT_by_tair.pdf", width = 3.45, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_minmax); dev.off()

# --- right: per-plant slope of ΔT against air temperature --------------------
# one regression per plant across the three days
beta <- daily %>%
  group_by(id, cluster) %>%
  filter(n() >= 3, sd(Tair) > 0) %>%
  summarise(beta_min = coef(lm(deltaT_min ~ Tair))[2],
            beta_max = coef(lm(deltaT_max ~ Tair))[2], .groups = "drop")

beta_long <- bind_rows(
  beta %>% transmute(cluster, value = beta_min, stat = "min"),
  beta %>% transmute(cluster, value = beta_max, stat = "max")
)

beta_letters <- bind_rows(lapply(c("min", "max"), function(st) {
  sub <- beta_long[beta_long$stat == st, ]
  fit <- aov(value ~ cluster, data = sub)
  pv  <- TukeyHSD(fit)[[1]][, "p adj"]
  names(pv) <- rownames(TukeyHSD(fit)[[1]])
  L <- multcompLetters(pv)$Letters
  ypos <- tapply(sub$value, sub$cluster, quantile, 0.95, na.rm = TRUE)
  data.frame(cluster = names(L), Letters = unname(L), stat = st,
             y = ypos[names(L)] + if (st == "min") 0.25 else 0.2)
})) %>%
  mutate(cluster = factor(cluster, levels = names(pal)))

print(summary(aov(beta_max ~ cluster, data = beta)))
print(TukeyHSD(aov(beta_max ~ cluster, data = beta)))

p_beta <- ggplot(beta_long, aes(cluster, value, fill = cluster)) +
  geom_boxplot(data = subset(beta_long, stat == "min"), outlier.shape = NA,
               width = 0.5, color = "grey30", alpha = 0.6) +
  geom_boxplot(data = subset(beta_long, stat == "max"), outlier.shape = NA,
               width = 0.5, color = "grey30") +
  geom_text(data = beta_letters, aes(cluster, y, label = Letters),
            inherit.aes = FALSE, size = 5, fontface = "bold",
            family = "Helvetica") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x)) +
  facet_wrap(~factor(stat, levels = c("min", "max")), strip.position = "bottom",
             labeller = as_labeller(c("min" = "Abs Min", "max" = "Abs Max"))) +
  labs(x = NULL, y = expression(beta ~ (Delta * italic(T) ~ "~" ~ italic(T)[air])),
       fill = "Cluster") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    strip.background  = element_blank(),
    strip.placement   = "outside",
    strip.text        = element_text(size = 14),
    legend.position   = "none",
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 10),
    axis.title        = element_text(size = 14),
    axis.title.y.left = element_text(margin = margin(r = -2)),
    panel.spacing     = unit(0.2, "cm"),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA)
  )

pdf("fig5a_right_beta.pdf", width = 2.75, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p_beta); dev.off()
