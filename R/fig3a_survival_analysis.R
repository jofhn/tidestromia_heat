# Figure 3A — survival by thermal cluster
# Input: data/fig3abcd_tidy_temps.csv (cluster assignments, fustat, time_to_four)

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

pal <- c("#E63946", "#F4A261", "#457B9D")

d <- read.csv("data/fig3abcd_tidy_temps.csv")[-1]

# k-means labels are arbitrary; reorder so cluster 1 is the fastest to die
d$cluster <- factor(d$cluster, levels = c("1", "3", "2"), labels = c("1", "2", "3"))

fit <- survfit(Surv(time_to_four, fustat) ~ cluster, data = d)
pv  <- surv_pvalue(fit, data = d)$pval.txt

surv <- ggsurvplot(
  fit, data = d,
  pval = FALSE, legend = "top", legend.title = "Cluster",
  legend.labs = c("1", "2", "3"),
  xlab = "Days to death", ylab = "Survival probability",
  palette = pal,
  font.family = "Helvetica",
  ggtheme = theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(
      legend.position   = "top",
      legend.title      = element_text(size = 13),
      legend.text       = element_text(size = 12),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.margin     = margin(b = -5),
      legend.box.margin = margin(b = -5),
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA)
    )
)

surv$plot <- surv$plot +
  annotate("text", x = 0, y = 0.15, label = pv, hjust = 0,
           size = 4, colour = "grey30", family = "Helvetica")

pdf("fig3a_survival.pdf", width = 3.5, height = 3,
    useDingbats = FALSE, bg = "transparent")
print(surv$plot)
dev.off()