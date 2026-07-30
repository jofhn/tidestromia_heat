# Figure 2B — QQ plot of the DV/Perimeter association test
# Inputs: data/DVP_imputed_sorted_100kb_speciesrange_ne.prune.prune.in
#         glm_tt4_imputed_dvp_allchr_pc1-25_100kb_speciesrange_ne.tsv (1.4 GB, not deposited; regenerate per wgs/survival_panel.md, section 13)
#
# Association results are from PLINK 2 (see wgs/survival_panel.md, section 13).

library(data.table)
library(ggplot2)

# --- data --------------------------------------------------------------------
glm <- fread("/mnt/gs21/scratch/feehanj1/analysis/25363Srh/imputation/glm/dvp/glm_tt4_imputed_dvp_allchr_pc1-25_100kb_speciesrange_ne.tsv",
             select = "P")

# variants with NA P failed the VIF check and are dropped
p <- glm$P
p <- p[is.finite(p) & p > 0 & p <= 1]
n <- length(p)

# --- genomic inflation and significance threshold ----------------------------
chisq     <- qchisq(1 - p, df = 1)
lambda_gc <- median(chisq) / qchisq(0.5, df = 1)

# Bonferroni over LD-pruned (approximately independent) variants
n_independent <- nrow(fread("data/DVP_imputed_sorted_100kb_speciesrange_ne.prune.prune.in", header = FALSE))
bonferroni <- 0.05 / n_independent

message("variants tested: ", n)
message("independent variants: ", n_independent)
message("lambda_GC: ", round(lambda_gc, 3))
message("Bonferroni threshold: ", signif(bonferroni, 3))

# --- plot --------------------------------------------------------------------
qq_df <- data.frame(
  exp      = -log10(ppoints(n)),
  obs      = -log10(sort(p)),
  ci_upper = -log10(qbeta(0.025, 1:n, n - 1:n + 1)),
  ci_lower = -log10(qbeta(0.975, 1:n, n - 1:n + 1))
)

# keep all points above -log10(P) = 2, subsample the dense null region
keep <- qq_df$obs > 2 | seq_len(n) %% 100 == 0
qq_df <- qq_df[keep, ]

qq <- ggplot(qq_df, aes(exp, obs)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "grey50", alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 0.7, shape = 16) +
  labs(x = expression(Expected ~ -log[10](italic(P))),
       y = expression(Observed ~ -log[10](italic(P)))) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    axis.text        = element_text(size = 10),
    axis.title       = element_text(size = 14),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("fig2b_qq.pdf", width = 2.65, height = 2.75, useDingbats = FALSE)
print(qq)
dev.off()
