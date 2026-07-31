# Figure 2C — Manhattan plot of the DV/Perimeter association test
# Inputs: data/DVP_imputed_sorted_100kb_speciesrange_ne.prune.prune.in
#         data/ld_blocks_full.tsv (LD block boundaries, see fig2d)
#         glm_dvp_full.tsv.gz  (1.4 GB uncompressed; in the Zenodo deposit; see section 13 in wgs/survival_panel.md)

library(data.table)
library(ggplot2)
library(ragg)

shift <- data.table::shift

# --- data --------------------------------------------------------------------
glm <- fread("data/glm_tt4_imputed_dvp_allchr_pc1-25_100kb_speciesrange_ne.tsv")
setnames(glm, "#CHROM", "CHROM")

n_independent <- nrow(fread("data/DVP_imputed_sorted_100kb_speciesrange_ne.prune.prune.in",
                            header = FALSE))
bnf <- -log10(0.05 / n_independent)

# --- cumulative position for the x-axis --------------------------------------
setkey(glm, CHROM, POS)
glm[, chr_n := as.integer(sub("Chr", "", CHROM))]

chr_offsets <- glm[, .(chr_len = max(POS)), by = chr_n][order(chr_n)]
chr_offsets[, offset := cumsum(shift(chr_len, fill = 0))]

glm <- merge(glm, chr_offsets[, .(chr_n, offset)], by = "chr_n")
glm[, pos_cum := POS + offset]

axis_df <- glm[, .(center = mean(pos_cum)), by = chr_n][order(chr_n)]

# --- flag variants inside an LD block ----------------------------------------
ld_blocks <- fread("data/ld_blocks_full.tsv")
ld_blocks[, chr_n := as.integer(sub("Chr", "", chrom))]
ld_blocks <- merge(ld_blocks, chr_offsets[, .(chr_n, offset)], by = "chr_n")
ld_blocks[, `:=`(xmin = block_start + offset, xmax = block_end + offset)]
ld_blocks <- ld_blocks[pop == "DVP"]

# the Chr2 block is narrow; widened for visibility only
ld_blocks[label == "Chr2_178968191", `:=`(xmin = xmin - 5e5, xmax = xmax + 5e5)]

glm[, in_block := FALSE]
for (i in seq_len(nrow(ld_blocks))) {
  glm[pos_cum >= ld_blocks$xmin[i] & pos_cum <= ld_blocks$xmax[i],
      in_block := TRUE]
}

# --- plot --------------------------------------------------------------------
manhattan <- ggplot(glm[!is.na(P)], aes(pos_cum, -log10(P))) +
  geom_point(data = glm[!is.na(P) & !in_block],
             aes(color = factor(chr_n %% 2)), size = 0.6) +
  geom_point(data = glm[!is.na(P) & in_block],
             color = "pink3", size = 1.1) +
  scale_color_manual(values = c("0" = "#7a8a78", "1" = "#a8b2a6"), guide = "none") +
  scale_x_continuous(breaks = axis_df$center,
                     labels = paste0("Chr", axis_df$chr_n)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  geom_hline(yintercept = bnf, linetype = "dotted", color = "red3") +
  labs(x = "Genomic position", y = expression(-log[10](italic(P)))) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 12),
    axis.title       = element_text(size = 14),
    axis.title.y     = element_text(size = 14, margin = margin(r = -2)),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

agg_png("fig2c_manhattan.png", width = 3.7, height = 2.75, units = "in", res = 300, background = "transparent")
print(manhattan)
dev.off()