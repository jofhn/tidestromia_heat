# DV/Perimeter locus boundaries with candidate genes (Chr1, Chr2, Chr6)
# Inputs: data/DVP_Chr1_84178924_ld.vcor
#         data/DVP_Chr2_178968191_ld.vcor
#         data/DVP_Chr6_94545774_ld_wide.vcor
#         data/candidate_genes.tsv
#
# LD was computed in PLINK 2 around each clumped index SNP.

library(data.table)
library(ggplot2)

# --- locus definitions -------------------------------------------------------
# r2 threshold and plotting buffer were set per locus by inspection
loci <- data.table(
  chrom     = c(1, 2, 6),
  pos       = c(84178924, 178968191, 94545774),
  file      = c("DVP_Chr1_84178924_ld.vcor",
                "DVP_Chr2_178968191_ld.vcor",
                "DVP_Chr6_94545774_ld_wide.vcor"),
  r2_thresh = c(0.4, 0.175, 0.1),
  buffer_bp = c(500e3, 100e3, 5000e3)
)

# --- LD block boundaries -----------------------------------------------------
# walk outward from the index SNP; the block ends where n_consecutive
# successive windows contain no variant above the r2 threshold
find_ld_boundary <- function(ld, index_pos, r2_threshold, window_kb = 100, n_consecutive = 3)
  {
  window_bp <- window_kb * 1000
  
  edge <- function(snps, forward) {
    boundary <- index_pos
    consecutive <- 0
    for (j in seq_len(nrow(snps))) {
      pos <- snps[j, POS_B]
      win <- if (forward) snps[POS_B >= pos & POS_B <= pos + window_bp]
      else         snps[POS_B >= pos - window_bp & POS_B <= pos]
      if (all(win$UNPHASED_R2 < r2_threshold)) {
        consecutive <- consecutive + 1
        if (consecutive >= n_consecutive) break
      } else {
        consecutive <- 0
        boundary <- pos
      }
    }
    boundary
  }
  
  list(start = edge(ld[POS_B < index_pos][order(-POS_B)], FALSE),
       end   = edge(ld[POS_B > index_pos][order(POS_B)],  TRUE))
}

# --- candidate genes ---------------------------------------------------------
# one row per gene: the best-supported ortholog assignment
cands <- fread("data/candidate_genes.tsv")
cands[, c("chr", "pos") := tstrsplit(index_snp, ":", keep = 1:2,
                                     type.convert = TRUE)]
gene_coords <- cands[, .SD[order(-as.integer(is_rbh == "TRUE"), -pident,
                                 -qcov, -bitscore, na.last = TRUE)][1],
                     by = tio_gene][, .(tio_gene, chr, gene_start, gene_end)]

# --- plot --------------------------------------------------------------------
plots <- lapply(seq_len(nrow(loci)), function(i) {
  ld    <- fread(file.path("data", loci$file[i]))
  block <- find_ld_boundary(ld, loci$pos[i], loci$r2_thresh[i])
  
  plot_min <- block$start - loci$buffer_bp[i]
  plot_max <- block$end   + loci$buffer_bp[i]
  
  ld_zoom <- ld[POS_B >= plot_min & POS_B <= plot_max]
  genes   <- gene_coords[chr == loci$chrom[i] &
                           gene_end >= plot_min & gene_start <= plot_max]
  
  ggplot(ld_zoom, aes(POS_B / 1e6, UNPHASED_R2, color = UNPHASED_R2)) +
    geom_rect(data = genes,
              aes(xmin = gene_start / 1e6, xmax = gene_end / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "#E9C46A", alpha = 0.8, inherit.aes = FALSE) +
    geom_point(size = 2, alpha = 0.5) +
    geom_vline(xintercept = loci$pos[i] / 1e6, color = "red3",
               linetype = "dotted") +
    geom_vline(xintercept = c(block$start, block$end) / 1e6,
               color = "grey20", linetype = "dotted") +
    coord_cartesian(xlim = c(plot_min, plot_max) / 1e6, ylim = c(0, 1)) +
    scale_color_gradientn(colors = c("grey80", "#A5C8D0", "#1F4E5F"),
                          limits = c(0, 1), oob = scales::squish,
                          na.value = "grey80") +
    labs(x = paste0("Chromosome ", loci$chrom[i], " (Mb)"),
         y = expression(italic(r)^2)) +
    theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(size = 12, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 14),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    )
})

for (i in seq_len(nrow(loci))) {
  pdf(sprintf("fig2d_chr%d.pdf", loci$chrom[i]),
      width = 4, height = 2.5, useDingbats = FALSE, bg = "transparent")
  print(plots[[i]])
  dev.off()
}
