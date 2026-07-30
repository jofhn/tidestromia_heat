# Figure 1G — survival health status scores through the lethal heat trial
# Input: data/fig1g_survival_health_status_scores.csv (row, col, then one column per day d7-d29; scores 0-4, "-" = missing)

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

score_labels <- c("0" = "healthy", "1" = "burnt", "2" = "sick",
                  "3" = "very_sick", "4" = "dead")

day_levels <- paste0("d", 8:29)

# --- data --------------------------------------------------------------------
# Score 0-4 per plant per day; missing values excluded. Counts are normalised
# within each day so bars sum to 1. Day 7 is dropped (pre-treatment).
combined <- read.csv("data/fig1g_survival_health_status_scores.csv") %>%
  pivot_longer(-c(row, col), names_to = "day", values_to = "score",
               values_transform = list(score = as.character)) %>%
  filter(score %in% names(score_labels)) %>%
  mutate(survivability = score_labels[score]) %>%
  count(day, survivability, name = "sums") %>%
  filter(day %in% day_levels) %>%
  group_by(day) %>%
  mutate(normsums = sums / sum(sums)) %>%
  ungroup() %>%
  mutate(
    day = factor(day, levels = day_levels),
    survivability = factor(survivability, levels = unname(score_labels))
  )

# --- plot --------------------------------------------------------------------
plot <- ggplot(combined, aes(x = day, y = normsums, fill = survivability)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2, width = 0.95) +
  scale_fill_manual(
    values = c(healthy = "#4A6741", burnt = "#7A9668", sick = "#A8BF96",
               very_sick = "#DEB8C0", dead = "#C47888"),
    labels = c(healthy = "Alive 1", burnt = "Alive 2", sick = "Alive 3",
               very_sick = "Dying", dead = "Dead")
  ) +
  scale_x_discrete(labels = gsub("d", "", day_levels)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
  labs(x = "Days Post-Germination", y = NULL, fill = "Score") +
  # axes drawn manually so the treatment timeline can sit above the panel
  annotate("segment", x = 0.5, xend = 0.5,  y = 0, yend = 1, linewidth = 0.5) +
  annotate("segment", x = 0.5, xend = 22.5, y = 0, yend = 0, linewidth = 0.5) +
  # duration bars (hours at each stage)
  annotate("segment", x = c(0.75, 7.75, 17.75), xend = c(7.25, 17.25, 22.25),
           y = 1.04, yend = 1.04, color = "grey40", linewidth = 0.6) +
  annotate("text", x = c(4, 12.5, 20, 22.6), y = 1.15,
           label = c("6", "8", "24", "hr"), size = 4.5, family = "Helvetica") +
  # temperature row
  annotate("segment", x = c(4.75, 17.75), xend = c(17.25, 22.25),
           y = 1.26, yend = 1.26, color = "grey40", linewidth = 0.6) +
  annotate("text", x = 1:4, y = 1.36, label = c("37", "44", "51", "55"),
           size = 3.5, family = "Helvetica") +
  annotate("text", x = c(11, 20, 22.6), y = 1.36, label = c("60", "30", "°C"),
           size = 4.5, family = "Helvetica") +
  # phase row
  annotate("segment", x = c(0.75, 4.75), xend = c(4.25, 17.25),
           y = 1.47, yend = 1.47, color = "grey40", linewidth = 0.6) +
  annotate("text", x = c(2.5, 11), y = 1.59,
           label = c("Acclimate", "Heat Treatment"), size = 5,
           family = "Helvetica") +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position  = "none",
    axis.line        = element_blank(),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 10, angle = 30),
    axis.title.x     = element_text(size = 14),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    plot.margin      = margin(t = 10, r = 10, b = 2, l = 20)
  )

# y-axis title placed manually, since the panel extends above y = 1
draw_plot <- function() {
  print(plot)
  grid.text("Normalised Score", x = unit(0.025, "npc"), y = unit(0.45, "npc"),
            rot = 90, gp = gpar(fontfamily = "Helvetica", fontsize = 12))
}

pdf("fig1g_survival_health_status_scores.pdf", width = 4.75, height = 2.5, useDingbats = FALSE, bg = "transparent")
draw_plot()
dev.off()
