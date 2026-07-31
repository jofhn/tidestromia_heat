# Figure 5D — effect of ABA on survival of the painted leaf
# Inputs: data/fig5d_survival_scores.csv   (daily scores, treated and control leaf)
#         data/fig5c_treatments.csv        (ABA / mock / untreated)

library(dplyr)
library(ggplot2)
library(multcompView)
library(FSA)

PREDAYS <- 9    # last day at 30 °C; scoring days are counted from here

pal <- c(aba = "#E63946", mock = "#8B7BAD", untreated = "#A0785A")
lab <- c(aba = "ABA", mock = "Mock", untreated = "Untreated")

# --- data --------------------------------------------------------------------
# One leaf per plant was painted (columns dNT) and another left as an internal
# control (dNC). Scores run 0 (healthy) to 4 (dead).
scores <- read.csv("data/fig5d_survival_scores.csv")
treatments <- read.csv("data/fig5c_treatments.csv")

d <- merge(scores, treatments, by = c("tray", "row", "col", "cell", "id")) %>%
  rename(treatment = round) %>%
  filter(treatment %in% names(pal))     # drops the "fca" arm

# days to a score of 4; leaves still alive at the end are given the day after
# scoring stopped, so every leaf has an event
time_to_death <- function(df, suffix) {
  days <- grep(paste0("^d[0-9]+", suffix, "$"), names(df), value = TRUE)
  nums <- as.integer(sub(paste0("^d([0-9]+)", suffix, "$"), "\\1", days))
  days <- days[order(nums)]
  nums <- sort(nums)
  
  apply(as.matrix(df[, days]), 1, function(x) {
    i <- which(x == 4)[1]
    (if (is.na(i)) max(nums) + 1 else nums[i]) - PREDAYS
  })
}

comp <- d %>%
  mutate(t4_treated = time_to_death(d, "T"),
         t4_control = time_to_death(d, "C"),
         delta      = t4_treated - t4_control,
         treatment  = factor(treatment, levels = names(pal))) %>%
  select(id, treatment, t4_treated, t4_control, delta)

# --- statistics --------------------------------------------------------------
# Residuals are not normal, so differences between treatments are tested by
# Kruskal-Wallis with Dunn's post-hoc rather than ANOVA.
print(shapiro.test(residuals(aov(delta ~ treatment, data = comp))))
print(kruskal.test(delta ~ treatment, data = comp))

dunn <- dunnTest(delta ~ treatment, data = comp, method = "bonferroni")$res
print(dunn)

pv <- setNames(dunn$P.adj, gsub(" - ", "-", dunn$Comparison))
letters_df <- data.frame(Letters = multcompLetters(pv)$Letters)
letters_df$treatment <- factor(rownames(letters_df), levels = names(pal))
letters_df$y <- max(comp$delta, na.rm = TRUE) + 1.6

# --- plot --------------------------------------------------------------------
p <- ggplot(comp, aes(treatment, delta, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "grey30", linewidth = 0.4) +
  geom_boxplot(width = 0.06, outlier.shape = NA, fill = "white",
               alpha = 0.4, color = "grey30", linewidth = 0.4) +
  geom_jitter(aes(color = treatment), width = 0.1, height = 0,
              size = 1.2, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_text(data = letters_df, aes(treatment, y, label = Letters),
            inherit.aes = FALSE, size = 5, fontface = "bold",
            family = "Helvetica") +
  scale_fill_manual(values = pal, labels = lab) +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = lab) +
  scale_y_continuous(breaks = seq(-13, 5, 2),
                     labels = function(x) gsub("-", "\u2212", x)) +
  labs(x = "Treatment", y = expression(Delta * " Survival (Days)")) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position   = "none",
    axis.text.x       = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y       = element_text(size = 12),
    axis.title        = element_text(size = 14),
    axis.title.y.left = element_text(margin = margin(r = -5)),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA)
  )

pdf("fig5d_aba_survival.pdf", width = 2.15, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(p); dev.off()