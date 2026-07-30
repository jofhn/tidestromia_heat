# Figure 3D — predictors of death at 60 °C
# Inputs: data/fig3d_tidy_counts.csv                     (scores, survival status)
#         data/fig3bcd_radiuses_leaves_combined_all.csv    (IR camera leaf temperatures)
#         data/fig3d_leaf_measurements_d16.csv           (leaf area, shape)
#         data/fig3d_leaf_greenness.csv                  (green ratio)

library(dplyr)
library(ggplot2)
library(MCMCglmm)

TARGET_DATE <- as.POSIXct("2025-04-16", tz = "UTC")   # 60 °C day
WIN_LO <- 1030
WIN_HI <- 1500

# --- plants ------------------------------------------------------------------
# Exclude cell 14 (a small cluster), plants scored 5 (sick) on d15-d16, and
# plants scored 4 (unestablished) before the heat treatment on d12-d15.
scores <- read.csv("data/fig3d_tidy_counts.csv")[-1] %>%
  filter(cell != 14) %>%
  filter(if_all(c(d15, d16), ~ .x != "5")) %>%
  filter(if_all(c(d12, d13, d14, d15), ~ .x != "4")) %>%
  select(id, tray, row, col, cell, individual, fustat)

# --- leaf traits -------------------------------------------------------------
measurements <- read.csv("data/fig3d_leaf_measurements_d16.csv") %>%
  rename(cell = 1, area = 2, perimeter = 3, length = 4, width = 5) %>%
  filter(area > 0, area <= 10) %>%
  select(cell, area, perimeter, length, width, AR)

greenness <- read.csv("data/fig3d_leaf_greenness.csv")

# --- leaf temperature, averaged over the analysis window ---------------------
ir <- read.csv("data/fig3bcd_radiuses_leaves_combined_all.csv")
dt <- as.POSIXct(gsub("Record_", "", ir$csv), format = "%Y-%m-%d_%H-%M-%S", tz = "UTC")
ir$csv_date     <- as.POSIXct(format(dt, "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC")
ir$numeric_time <- as.numeric(format(dt, "%H%M"))

leaf_temp <- ir %>%
  filter(csv_date == TARGET_DATE, numeric_time >= WIN_LO, numeric_time <= WIN_HI) %>%
  group_by(cell) %>%
  summarise(avg_leaf_temp = mean(area_mean_3, na.rm = TRUE), .groups = "drop")

# --- assemble ----------------------------------------------------------------
d <- scores %>%
  inner_join(measurements, by = "cell") %>%
  inner_join(greenness,    by = "cell") %>%
  inner_join(leaf_temp,    by = "cell")

# distance to the nearest edge of the 9 x 18 tray
d$dist_edge <- apply(d[, c("row", "col")], 1,
                     function(x) min(x[1] - 1, 9 - x[1], x[2] - 1, 18 - x[2]))

# number of occupied neighbouring cells in the same tray (king's move)
d$local_density <- sapply(seq_len(nrow(d)), function(i) {
  sum(d$tray == d$tray[i] &
        abs(d$row - d$row[i]) <= 1 &
        abs(d$col - d$col[i]) <= 1) - 1
})

# predictors standardised so effects are comparable
d <- d %>%
  mutate(across(c(avg_leaf_temp, area, AR, green_ratio, local_density, dist_edge),
                ~ as.numeric(scale(.x)),
                .names = "{.col}_z"),
         tray       = as.factor(tray),
         individual = as.factor(individual))

# --- model -------------------------------------------------------------------
set.seed(42)
model <- MCMCglmm(
  fixed  = fustat ~ area_z + AR_z + avg_leaf_temp_z + green_ratio_z +
    local_density_z + dist_edge_z,
  random = ~ tray + individual,
  family = "categorical",
  data   = d,
  prior  = list(R = list(V = 1, fix = 1),
                G = list(G1 = list(V = 1, nu = 1),
                         G2 = list(V = 1, nu = 1))),
  pr = FALSE, nitt = 100000, burnin = 10000, thin = 10, verbose = FALSE
)
print(summary(model))

# convergence diagnostics
print(autocorr.diag(model$Sol))
print(effectiveSize(model$Sol))
print(heidel.diag(model$Sol))

# --- plot --------------------------------------------------------------------
name_map <- c(area_z            = "Leaf Area",
              AR_z              = "Leaf Shape",
              avg_leaf_temp_z   = "Leaf Temperature",
              green_ratio_z     = "Green Ratio",
              local_density_z   = "Local Density",
              dist_edge_z       = "Distance to Edge")

s <- summary(model)$solutions
est <- data.frame(
  term      = rownames(s),
  post.mean = s[, "post.mean"],
  l95       = s[, "l-95% CI"],
  u95       = s[, "u-95% CI"],
  pMCMC     = s[, "pMCMC"]
) %>%
  filter(term != "(Intercept)") %>%
  mutate(signif = cut(pMCMC, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                      labels = c("***", "**", "*", ".", "")),
         term   = name_map[term])

p <- ggplot(est, aes(post.mean, reorder(term, post.mean))) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_errorbarh(aes(xmin = l95, xmax = u95), width = 0.2,
                 linewidth = 0.5, color = "grey30") +
  geom_point(size = 2) +
  geom_text(aes(label = signif), vjust = -0.05, size = 5, family = "Helvetica") +
  labs(x = "Effect on log-odds of death", y = "Predictors") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    axis.text        = element_text(size = 10),
    axis.title.x     = element_text(size = 12),
    axis.title.y     = element_text(size = 14, vjust = -1.5),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("fig3d_MCMCglmm.pdf", width = 3.75, height = 2.5,
    useDingbats = FALSE, bg = "transparent")
print(p)
dev.off()
