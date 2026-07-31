# Figure 4 and figs. S12-S13 — RNA-seq of the heat series
# Inputs: data/fig4_star_gene_counts_unstranded.tsv   (gene x sample counts)
#         data/fig4_RNAseq_hotcold_plan.csv           (sample metadata)
#         data/fig4_gene_families_annotation.xlsx     (GO terms, Arabidopsis orthologs)
#         data/dvp_candidates.tsv                     (GWAS candidate genes per locus)
#
# Counts are from STAR --quantMode GeneCounts (see rnaseq/hotcold.md).

library(DESeq2)
library(apeglm)
library(ashr)
library(RUVSeq)
library(data.table)
library(readxl)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(eulerr)

LEVS      <- c("30","37","44","51","55_cold","55_hot",
               "d9","d10","d11","d12","d13")
HEAT_LEVS <- c("37","44","51","55_cold","55_hot")
DEV_MATCH <- c("37"="d10", "44"="d11", "51"="d12",
               "55_cold"="d13", "55_hot"="d13")
LFC_THR   <- log2(1.5)

col_map <- c("30"="#A7C957","37"="#6A994E","44"="#2EC4B6","51"="#F4A261",
             "55_cold"="#457B9D","55_hot"="#E63946",
             "d9"="#000000","d10"="#343A40","d11"="#ADB5BD",
             "d12"="#CED4DA","d13"="#E9ECEF")

col_37 <- "#F1C0A0"; col_44 <- "#EAA25E"; col_51 <- "#F4A261"
col_cold <- "#457B9D"; col_hot <- "#E63946"; col_q2a <- "grey25"

# --- counts and metadata -----------------------------------------------------
counts <- read.delim("data/fig4_star_gene_counts_unstranded.tsv",
                     check.names = FALSE, row.names = 1)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# unplaced scaffolds are excluded
counts <- counts[!grepl("^TIOBLUn", rownames(counts)), ]

# JF11 and JF18 were mislabelled at library prep; the swap was identified by
# leave-one-out centroid correlation on VST counts and confirmed by PCA
swap <- match(c("JF11", "JF18"), colnames(counts))
colnames(counts)[swap] <- colnames(counts)[rev(swap)]

samples <- read.csv("data/fig4_RNAseq_hotcold_plan.csv",
                    stringsAsFactors = FALSE)
samples$novogene <- paste0("JF", sprintf("%02d", as.integer(samples$JF)))
samples <- samples[samples$novogene %in% colnames(counts), ]
rownames(samples) <- samples$novogene
samples <- samples[colnames(counts), ]
samples$sample <- factor(samples$sample, levels = LEVS)

# --- panel 4A: PCA on all samples --------------------------------------------
vsd <- vst(DESeqDataSetFromMatrix(counts, samples, ~ 1), blind = TRUE)
pca <- prcomp(t(assay(vsd)))
pct <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     sample = samples$sample)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = sample)) +
  geom_point(size = 4, alpha = 0.6) +
  scale_color_manual(values = col_map, breaks = LEVS) +
  labs(x = paste0("PC1 (", round(pct[1], 1), "%)"),
       y = paste0("PC2 (", round(pct[2], 1), "%)"),
       color = "Timepoint") +
  guides(color = guide_legend(ncol = 2, title.hjust = 0.45)) +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.key.size   = unit(0.3, "cm"),
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 13),
    axis.text         = element_text(size = 10),
    axis.title        = element_text(size = 14),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA)
  )

pdf("fig4a_pca.pdf", width = 4.4, height = 2.2,
    useDingbats = FALSE, bg = "transparent")
print(p_pca); dev.off()

# --- differential expression -------------------------------------------------
# keep genes with >=10 reads in >=3 samples (the smallest group size)
counts <- counts[rowSums(counts >= 10) >= 3, ]

dds <- DESeq(DESeqDataSetFromMatrix(counts, samples, ~ sample))

# contrast vector against resultsNames(); intercept is the "30" level
build_contrast <- function(dds, coefs) {
  v <- setNames(rep(0, length(resultsNames(dds))), resultsNames(dds))
  for (nm in names(coefs)) v[paste0("sample_", nm, "_vs_30")] <- coefs[[nm]]
  v
}

# Q1: each heat timepoint vs its developmental control, relative to baseline
q1 <- rbindlist(lapply(HEAT_LEVS, function(tp) {
  v <- build_contrast(dds, setNames(list(1, -1, 1),
                                    c(tp, DEV_MATCH[[tp]], "d9")))
  thr <- results(dds, contrast = v, lfcThreshold = LFC_THR,
                 altHypothesis = "greaterAbs", alpha = 0.05)
  shr <- lfcShrink(dds, contrast = v, type = "ashr", lfcThreshold = LFC_THR)
  data.table(gene = rownames(thr),
             log2FoldChange = shr$log2FoldChange,
             padj = thr$padj,
             contrast = paste0("(", tp, "-", DEV_MATCH[[tp]], ")-(30-d9)"))
}))
q1[, sig := !is.na(padj) & padj < 0.05]

# Q2a: 55_cold vs 55_hot directly, releveled so apeglm can shrink the coefficient
dds_55 <- dds
dds_55$sample <- relevel(dds_55$sample, ref = "55_hot")
dds_55 <- nbinomWaldTest(dds_55)

q2a <- as.data.table(as.data.frame(
  lfcShrink(dds_55, coef = "sample_55_cold_vs_55_hot",
            type = "apeglm", lfcThreshold = LFC_THR)
), keep.rownames = "gene")
q2a[, sig := svalue < 0.005 & abs(log2FoldChange) > LFC_THR]

print(q1[, .(n_sig = sum(sig)), by = contrast])
cat("Q2a sig:", sum(q2a$sig, na.rm = TRUE), "\n")

# --- RUVs check --------------------------------------------------------------
# Refit with two factors of unwanted variation; reported in the methods as a
# consistency check, not used for any figure.
groups <- split(seq_len(ncol(counts)), samples$sample)
scIdx  <- do.call(rbind, lapply(groups, function(x)
  c(x, rep(-1, max(lengths(groups)) - length(x)))))

set <- RUVs(newSeqExpressionSet(counts,
                                phenoData = data.frame(samples,
                                                       row.names = colnames(counts))),
            cIdx = rownames(counts), k = 2, scIdx = scIdx)

meta_ruv <- cbind(samples, pData(set)[, c("W_1", "W_2")])
dds_ruv  <- DESeq(DESeqDataSetFromMatrix(counts, meta_ruv,
                                         ~ W_1 + W_2 + sample))

# W factors should not track biological condition
mm_pc1 <- prcomp(model.matrix(~ sample, data = meta_ruv)[, -1], scale. = TRUE)$x[, 1]
cat("W_1 vs condition PC1: r =", round(cor(meta_ruv$W_1, mm_pc1), 3), "\n")
cat("W_2 vs condition PC1: r =", round(cor(meta_ruv$W_2, mm_pc1), 3), "\n")

q1_ruv <- rbindlist(lapply(HEAT_LEVS, function(tp) {
  v <- build_contrast(dds_ruv, setNames(list(1, -1, 1),
                                        c(tp, DEV_MATCH[[tp]], "d9")))
  res <- results(dds_ruv, contrast = v, alpha = 0.05)
  data.table(gene = rownames(res), log2FoldChange = res$log2FoldChange,
             pvalue = res$pvalue,
             contrast = paste0("(", tp, "-", DEV_MATCH[[tp]], ")-(30-d9)"))
}))
q1_ruv[, q := p.adjust(pvalue, method = "BH"), by = contrast]
q1_ruv[, sig := !is.na(q) & q < 0.05]

for (cn in unique(q1$contrast)) {
  a <- q1[contrast == cn & sig, gene]; b <- q1_ruv[contrast == cn & sig, gene]
  cat(sprintf("  %-25s std=%5d ruv=%5d jaccard=%.3f\n", cn, length(a), length(b),
              length(intersect(a, b)) / length(union(a, b))))
}

# --- developmentally corrected fold changes, per timepoint -------------------
log_mat <- log2(counts(dds, normalized = TRUE) + 1)
cond_mean <- function(cond) rowMeans(log_mat[, samples$sample == cond, drop = FALSE])

baseline  <- cond_mean("30") - cond_mean("d9")
resid_mat <- sapply(HEAT_LEVS, function(tp)
  (cond_mean(tp) - cond_mean(DEV_MATCH[[tp]])) - baseline)

# --- DE gene sets ------------------------------------------------------------
sig_by_tp <- function(tp) q1[contrast == paste0("(", tp, "-", DEV_MATCH[[tp]],
                                                ")-(30-d9)") & sig, gene]

sig_37 <- sig_by_tp("37"); sig_44 <- sig_by_tp("44"); sig_51 <- sig_by_tp("51")
sig_cold <- sig_by_tp("55_cold"); sig_hot <- sig_by_tp("55_hot")

up_hot  <- q2a[sig & log2FoldChange < 0, gene]
up_cold <- q2a[sig & log2FoldChange > 0, gene]
sig_q2a <- q2a[sig == TRUE, gene]

de_sets <- list(
  "37"                = sig_37,
  "44"                = sig_44,
  "51"                = sig_51,
  "55_cold"           = sig_cold,
  "55_cold_unique"    = setdiff(sig_cold, sig_hot),
  "55_hot"            = sig_hot,
  "55_hot_unique"     = setdiff(sig_hot, sig_cold),
  "55_cold vs 55_hot" = sig_q2a
)

# --- panel 4B: heatmap of the 55 °C contrast ---------------------------------
set.seed(42)
bg <- sample(intersect(q2a[sig == FALSE, gene], rownames(resid_mat)), 600)

slices <- factor(c(rep("Up in 55_hot",      length(up_hot)),
                   rep("Up in 55_cold",     length(up_cold)),
                   rep("Background (n.s.)", length(bg))),
                 levels = c("Up in 55_hot", "Up in 55_cold", "Background (n.s.)"))

fill_lim <- quantile(abs(resid_mat), 0.99, na.rm = TRUE)
col_fun  <- colorRamp2(c(-fill_lim, 0, fill_lim),
                       c("#2C5F8D", "#FAF7F2", "#D4A017"))

ht <- Heatmap(
  resid_mat[c(up_hot, up_cold, bg), ],
  name = "log2FC", col = col_fun,
  cluster_rows = TRUE, clustering_method_rows = "ward.D2",
  cluster_columns = FALSE, cluster_row_slices = FALSE,
  row_split = slices, row_gap = unit(2, "mm"),
  show_row_names = FALSE, column_names_rot = 30,
  column_names_gp = gpar(fontsize = 13, fontfamily = "Helvetica"),
  column_title = "Timepoint", column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontfamily = "Helvetica"),
  row_title_gp = gpar(fontsize = 14, fontfamily = "Helvetica"),
  row_title_rot = 90,
  heatmap_legend_param = list(
    title = expression(log[2] * " FC"),
    title_gp = gpar(fontsize = 12, fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
    title_position = "topcenter", direction = "horizontal",
    legend_width = unit(4, "cm")),
  rect_gp = gpar(col = NA), use_raster = FALSE
)

pdf("fig4b_heatmap.pdf", width = 4, height = 9,
    useDingbats = FALSE, bg = "transparent")
draw(ht, heatmap_legend_side = "top", background = "transparent")
dev.off()

# --- GO annotation -----------------------------------------------------------
# GO terms and Arabidopsis orthologs both come from the gene-family table
matts <- as.data.table(read_excel("data/fig4_gene_families_annotation.xlsx"))
setnames(matts, gsub("\r\n", " ", colnames(matts)))

ROOTS <- c("GO:0008150", "GO:0003674", "GO:0005575")

parse_go <- function(col, label) {
  d <- matts[get(col) != "NA", .(GeneID, raw = get(col))]
  d[, .(go_id = unlist(str_extract_all(raw, "GO:\\d{7}"))), by = GeneID][
    , ontology := label][]
}

go_long <- rbindlist(list(
  parse_go("GO Terms (Biological Process)",   "BP"),
  parse_go("GO Terms (Molecular Function)",   "MF"),
  parse_go("GO Terms (Cellular Compartment)", "CC")
))
go_long <- unique(go_long[!go_id %in% ROOTS], by = c("GeneID", "go_id"))

term2name <- rbindlist(lapply(
  c("GO Terms (Biological Process)", "GO Terms (Molecular Function)",
    "GO Terms (Cellular Compartment)"),
  function(col) matts[get(col) != "NA", .(raw = get(col))]))
term2name <- term2name[, .(entry = unlist(strsplit(raw, "; "))), by = .I][
  , .(go_id   = str_extract(entry, "GO:\\d{7}"),
      go_name = str_extract(entry, "(?<=\\{).*(?=\\})"))]
term2name <- unique(term2name[!is.na(go_id) & !go_id %in% ROOTS])

# background: expressed genes carrying any annotation
bg_genes <- intersect(rownames(counts), unique(go_long$GeneID))

run_enrich <- function(fg, ont, label) {
  res <- enricher(intersect(fg, bg_genes), universe = bg_genes,
                  TERM2GENE = go_long[ontology == ont, .(go_id, GeneID)],
                  TERM2NAME = term2name,
                  pvalueCutoff = 0.05, qvalueCutoff = 0.10,
                  minGSSize = 5, maxGSSize = 500)
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  as.data.table(res@result)[, `:=`(ontology = ont, gene_set = label)][]
}

fg_list <- list(Q2a_up_in_cold = up_cold, Q2a_up_in_hot = up_hot)

enrich_all <- rbindlist(lapply(names(fg_list), function(label) {
  rbindlist(lapply(c("BP", "MF", "CC"), run_enrich,
                   fg = fg_list[[label]], label = label), fill = TRUE)
}), fill = TRUE)

# --- panel 4C: GO dot plots --------------------------------------------------
parse_ratio <- function(x) sapply(strsplit(x, "/"),
                                  function(p) as.numeric(p[1]) / as.numeric(p[2]))

build_dot <- function(label) {
  dt <- enrich_all[gene_set == label & qvalue < 0.05 & Count >= 3,
                   .(ontology, Description, Count, qvalue,
                     fold_enrich = parse_ratio(GeneRatio) / parse_ratio(BgRatio))]
  dt[, Description := reorder(Description, -qvalue)][]
}

make_dot <- function(dt, title, x_max, cols) {
  ggplot(dt, aes(fold_enrich, Description, color = ontology)) +
    geom_point(aes(size = Count), alpha = 0.85) +
    geom_text(aes(label = Count), color = "white", size = 3, fontface = "bold") +
    facet_grid(ontology ~ ., scales = "free", space = "free_y") +
    scale_x_continuous(limits = c(0, x_max), oob = scales::squish) +
    scale_color_manual(values = cols, guide = "none") +
    scale_size_continuous(range = c(5, 9), guide = "none") +
    labs(x = "Fold Enrichment", y = NULL, title = title) +
    theme_classic(base_family = "Helvetica", base_size = 11) +
    theme(
      plot.title       = element_text(size = 16, hjust = 0.5),
      strip.background = element_blank(),
      strip.text       = element_text(size = 14),
      axis.text.y      = element_text(size = 12, lineheight = 0.8),
      axis.text.x      = element_text(size = 10),
      axis.title       = element_text(size = 14),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    )
}

p_hot <- make_dot(build_dot("Q2a_up_in_hot"), "Up in 55_Hot", 20,
                  c(BP = "#E8B7C5", MF = "#C56B8A", CC = "#7A2E3E"))
p_cold <- make_dot(build_dot("Q2a_up_in_cold"), "Up in 55_Cold", 30,
                   c(BP = "#A7C7A0", MF = "#5A9367", CC = "#2D4A2E"))

pdf("fig4c_go_hot.pdf",  width = 4.1, height = 5,   useDingbats = FALSE, bg = "transparent")
print(p_hot);  dev.off()
pdf("fig4c_go_cold.pdf", width = 4.1, height = 2.1, useDingbats = FALSE, bg = "transparent")
print(p_cold); dev.off()

# --- panels 4D and 4E: trajectories of enriched gene sets --------------------
plot_go_lines <- function(set_label, term, labels, title, cols, nudge) {
  genes <- enrich_all[gene_set == set_label & Description == term,
                      unlist(strsplit(geneID, "/"))]
  genes <- intersect(genes, rownames(resid_mat))
  stopifnot(length(genes) > 0)
  
  d <- melt(as.data.table(resid_mat[genes, , drop = FALSE], keep.rownames = "gene"),
            id.vars = "gene", variable.name = "timepoint", value.name = "residual")
  d[, timepoint := factor(timepoint, levels = HEAT_LEVS)]
  d[, short := labels[gene]]
  
  ggplot(d, aes(as.numeric(timepoint), residual, group = gene, color = gene)) +
    annotate("segment", x = 1, xend = 5, y = 0, yend = 0,
             linetype = "dotted", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_text_repel(data = d[timepoint == "55_hot" & !is.na(short)],
                    aes(label = short), direction = "y", hjust = 0,
                    nudge_x = nudge, xlim = c(5, NA), size = 4,
                    family = "Helvetica", segment.size = 0.3,
                    segment.color = "grey60", min.segment.length = 0,
                    show.legend = FALSE) +
    scale_x_continuous(breaks = seq_along(HEAT_LEVS), labels = HEAT_LEVS,
                       expand = expansion(mult = c(0.05, nudge / 1.5))) +
    annotate("segment", x = 1, xend = 5, y = -Inf, yend = -Inf, linewidth = 0.5) +
    scale_color_manual(values = colorRampPalette(cols)(length(genes))) +
    labs(x = "Timepoint", y = expression("log"[2] ~ "FC"), title = title) +
    theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(
      plot.title       = element_text(size = 14, hjust = 0.4),
      axis.text.x      = element_text(size = 12, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 12),
      axis.title.x     = element_text(size = 14),
      axis.title.y     = element_text(size = 14, margin = margin(r = -2)),
      axis.line.x      = element_blank(),
      legend.position  = "none",
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    )
}

p_photo <- plot_go_lines(
  "Q2a_up_in_hot", "photorespiration",
  c("TIOBL11G0064350" = "SHMT", "TIOBL05G0168430" = "HPR1"),
  "GO Photorespiration Genes",
  c("#E8B7C5", "#C56B8A", "#7A2E3E"), 0.25)

p_auxin <- plot_go_lines(
  "Q2a_up_in_cold", "auxin binding",
  c("TIOBL11G0074190" = "AUX/LAX", "TIOBL03G0121020" = "TIR/AFB",
    "TIOBL04G0142240" = "TIR/AFB", "TIOBL05G0184100" = "TIR/AFB"),
  "GO Auxin Binding Genes",
  c("#A7C7A0", "#5A9367", "#2D4226"), 0.5)

pdf("fig4d_photorespiration.pdf", width = 4.1, height = 3,
    useDingbats = FALSE, bg = "transparent")
print(p_photo); dev.off()

pdf("fig4e_auxin.pdf", width = 4.1, height = 3,
    useDingbats = FALSE, bg = "transparent")
print(p_auxin); dev.off()

# =============================================================================
# fig. S12
# =============================================================================

# --- S12B: pairwise overlap between DE gene sets -----------------------------
set_order <- names(de_sets)

label_cols <- c("37" = "#A7C957", "44" = "#6A994E", "51" = "#F4A261",
                "55_cold" = "#457B9D", "55_cold_unique" = "#2A4E66",
                "55_hot" = "#E63946", "55_hot_unique" = "#96202B",
                "55_cold vs 55_hot" = "grey20")

overlap_long <- rbindlist(lapply(set_order, function(s1) {
  rbindlist(lapply(set_order, function(s2) {
    data.table(set1 = s1, set2 = s2,
               count = length(intersect(de_sets[[s1]], de_sets[[s2]])))
  }))
}))
overlap_long[, denom := length(de_sets[[set1[1]]]), by = set1]
overlap_long[, fraction := round(count / denom, 2)]
overlap_long[, set1 := factor(set1, levels = set_order)]
overlap_long[, set2 := factor(set2, levels = set_order)]

p_matrix <- ggplot(overlap_long, aes(set2, set1, fill = fraction)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f\n(%d)", fraction, count)),
            color = ifelse(overlap_long$fraction > 0.5, "white", "grey20"),
            size = 2.1, lineheight = 0.85, family = "Helvetica") +
  scale_fill_gradient(low = "#FAF7F2", high = "#4A6741",
                      limits = c(0, 1), name = "Fraction") +
  scale_y_discrete(limits = rev(set_order)) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12,
                                      color = label_cols[set_order]),
    axis.text.y        = element_text(size = 12, color = label_cols[rev(set_order)]),
    axis.ticks         = element_blank(),
    axis.line          = element_blank(),
    legend.title       = element_text(size = 12, hjust = 0),
    legend.text        = element_text(size = 12, hjust = 0),
    legend.box.spacing = unit(0, "pt"),
    legend.key.width   = unit(0.3, "cm"),
    legend.key.height  = unit(0.5, "cm"),
    panel.background   = element_rect(fill = "transparent", color = NA),
    plot.background    = element_rect(fill = "transparent", color = NA),
    legend.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("figS12b_overlap_matrix.pdf", width = 5.25, height = 4.75, useDingbats = FALSE)
print(p_matrix); dev.off()

# --- S12C: DEG counts and Euler diagram --------------------------------------
sig_counts <- data.table(
  contrast = factor(c("37°C","44°C","51°C","55_cold","55_hot","55_cold vs 55_hot"),
                    levels = c("37°C","44°C","51°C","55_cold","55_hot",
                               "55_cold vs 55_hot")),
  n_sig = c(length(sig_37), length(sig_44), length(sig_51),
            length(sig_cold), length(sig_hot), length(sig_q2a))
)

p_bar <- ggplot(sig_counts, aes(contrast, n_sig, fill = contrast)) +
  geom_col() +
  geom_text(aes(label = n_sig), vjust = -0.3, size = 4, family = "Helvetica") +
  scale_fill_manual(values = c("37°C" = col_37, "44°C" = col_44, "51°C" = col_51,
                               "55_cold" = col_cold, "55_hot" = col_hot,
                               "55_cold vs 55_hot" = col_q2a)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Number of DEGs") +
  theme_classic(base_family = "Helvetica", base_size = 13) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 30, hjust = 1, size = 13),
    axis.text.y      = element_text(size = 11),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

pdf("figS12c_deg_counts.pdf", width = 2.85, height = 4.1,
    useDingbats = FALSE, bg = "transparent")
print(p_bar); dev.off()

fit_extreme <- euler(list("51°C" = sig_51, "55_cold" = sig_cold,
                          "55_hot" = sig_hot, "55_cold vs 55_hot" = sig_q2a),
                     shape = "ellipse")
cat("euler stress:", fit_extreme$stress,
    " diagError:", fit_extreme$diagError, "\n")

p_euler <- plot(fit_extreme,
                fills = list(fill = c(col_51, col_cold, col_hot, col_q2a),
                             alpha = 0.45),
                labels = FALSE,
                quantities = list(fontfamily = "Helvetica", cex = 0.7,
                                  col = "grey25"),
                edges = list(col = "grey40", lwd = 1),
                legend = list(fontfamily = "Helvetica", cex = 0.9, side = "right"))

pdf("figS12c_euler.pdf", width = 4, height = 4,
    useDingbats = FALSE, bg = "transparent")
print(p_euler); dev.off()

# --- S12D and S12E: GWAS candidate genes that are also DE --------------------
dvp_cands <- fread("data/candidate_genes.tsv")

LOCI <- c(Chr1 = "1:84178924:C:T",
          Chr2 = "2:178968191:A:G",
          Chr6 = "6:94545774:A:T")

gene_labels <- c(
  "TIOBL01G0015440" = "disordered;<br>transposase (Spm-like)",
  "TIOBL01G0015460" = "disordered;<br>Zn finger (CCHC)",
  "TIOBL01G0015560" = "PdBG4",
  "TIOBL01G0015600" = "GLX1",
  "TIOBL06G0198620" = "RING1B",
  "TIOBL06G0198630" = "TPX2",
  "TIOBL06G0198650" = "Reticulon family",
  "TIOBL06G0198660" = "PCMP-E75",
  "TIOBL06G0198690" = "unknown",
  "TIOBL06G0198700" = "RPN1A",
  "TIOBL06G0198740" = "proteophosphoglycan-<br>related",
  "TIOBL06G0198750" = "FLA family",
  "TIOBL06G0198910" = "AICARFT/IMPCHase",
  "TIOBL06G0198930" = "unknown",
  "TIOBL06G0199020" = "ECAP",
  "TIOBL06G0199260" = "disordered;<br>LEA domain",
  "TIOBL06G0199370" = "putative<br>F-box protein",
  "TIOBL06G0199390" = "unknown",
  "TIOBL06G0199410" = "aspartyl protease +<br>retrotransposon gag",
  "TIOBL06G0199550" = "Pentatricopeptide repeat"
)

# which candidate genes are DE, and in which contrasts
gwas_de <- rbindlist(c(
  lapply(HEAT_LEVS, function(tp)
    data.table(timepoint = tp,
               gene = intersect(sig_by_tp(tp), dvp_cands$tio_gene))),
  list(data.table(timepoint = "cold_v_hot",
                  gene = intersect(sig_q2a, dvp_cands$tio_gene)))
))
gwas_de <- merge(gwas_de, unique(dvp_cands[, .(gene = tio_gene, index_snp)]),
                 by = "gene")
gwas_de[, locus := names(LOCI)[match(index_snp, LOCI)]]

tp_cols <- c("37" = "#A7C957", "44" = "#6A994E", "51" = "#F4A261",
             "55_cold" = "#457B9D", "55_hot" = "#E63946",
             "cold_v_hot" = "grey20")

print(gwas_de[, .N, by = .(timepoint, locus)])

plot_locus <- function(chrom) {
  d <- gwas_de[locus == chrom & gene %in% rownames(resid_mat) &
                 !timepoint %in% c("37", "44")]
  if (nrow(d) == 0) return(invisible(NULL))
  d[, timepoint := factor(timepoint, levels = names(tp_cols))]
  d <- d[order(gene, timepoint)]
  
  # panel header: short name, gene ID, and the contrasts it is DE in
  d[, colored := paste0("<span style='color:", tp_cols[as.character(timepoint)],
                        ";'>", timepoint, "</span>")]
  labs <- d[, .(datasets = paste(colored, collapse = ", ")), by = gene]
  labs[, short := fifelse(is.na(gene_labels[gene]), gene, gene_labels[gene])]
  labs[, panel_label := paste0(short, "<br>", gene, "<br>DE: ", datasets)]
  
  rd <- rbindlist(lapply(unique(d$gene), function(g)
    data.table(gene = g, timepoint_x = factor(HEAT_LEVS, levels = HEAT_LEVS),
               residual = resid_mat[g, ])))
  rd <- merge(rd, labs[, .(gene, panel_label)], by = "gene")
  rd[, panel_label := factor(panel_label, levels = unique(panel_label))]
  
  ggplot(rd, aes(as.numeric(timepoint_x), residual, group = gene)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_line(linewidth = 0.8, color = "grey30") +
    geom_point(size = 1.8, color = "grey30") +
    scale_x_continuous(breaks = seq_along(HEAT_LEVS), labels = HEAT_LEVS) +
    facet_wrap(~ panel_label, ncol = 4) +
    labs(x = "Timepoint", y = expression("log"[2] ~ "FC"),
         title = bquote(.(chrom) ~ "genes" ~ intersect() ~ "DEGs")) +
    theme_classic(base_family = "Helvetica", base_size = 11) +
    theme(
      plot.title       = element_text(size = 13, hjust = 0),
      strip.background = element_blank(),
      strip.text       = element_markdown(size = 8, family = "Helvetica",
                                          lineheight = 1.2,
                                          margin = margin(t = 10, b = 0)),
      axis.text.x      = element_text(size = 10, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 12),
      panel.spacing    = unit(0.3, "lines"),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    )
}

pdf("figS12d_chr1_candidates.pdf", width = 5.75, height = 2.75,
    useDingbats = FALSE, bg = "transparent")
print(plot_locus("Chr1")); dev.off()

pdf("figS12e_chr6_candidates.pdf", width = 5.75, height = 7.5,
    useDingbats = FALSE, bg = "transparent")
print(plot_locus("Chr6")); dev.off()

# =============================================================================
# fig. S13 — core pathway homologs
# =============================================================================
# A TIOBL gene is called a homolog only where its FIRST Arabidopsis match is one
# of the canonical genes for that pathway.

core_family_plot <- function(canonical, at_to_family, fam_levels, go_term,
                             go_set, cols, title, ncol) {
  orth <- matts[sapply(`Arabidopsis Orthologs`,
                       function(x) any(grepl(paste(canonical, collapse = "|"), x))),
                .(gene = GeneID, ortholog = `Arabidopsis Orthologs`)]
  orth[, first_at := str_extract(ortholog, "AT\\dG\\d{5}")]
  hc <- orth[first_at %in% canonical]
  
  genes <- intersect(hc$gene, rownames(resid_mat))
  cat(title, "- high-confidence homologs:", length(genes), "\n")
  
  d <- melt(as.data.table(resid_mat[genes, , drop = FALSE], keep.rownames = "gene"),
            id.vars = "gene", variable.name = "timepoint", value.name = "residual")
  d[, timepoint := factor(timepoint, levels = HEAT_LEVS)]
  d <- merge(d, hc[, .(gene, first_at)], by = "gene")
  d[, family := factor(at_to_family[first_at], levels = fam_levels)]
  
  # genes recovered by the GO enrichment are highlighted
  in_go <- enrich_all[gene_set == go_set & Description == go_term,
                      unlist(strsplit(geneID, "/"))]
  d[, in_go := gene %in% in_go]
  
  ggplot(d, aes(as.numeric(timepoint), residual, group = gene, color = in_go)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq_along(HEAT_LEVS), labels = HEAT_LEVS) +
    scale_color_manual(values = cols) +
    facet_wrap(~ family, ncol = ncol) +
    labs(x = "Timepoint", y = expression("log"[2] ~ "FC"), title = title) +
    theme_classic(base_family = "Helvetica", base_size = 13) +
    theme(
      plot.title       = element_text(size = 14, hjust = 0),
      strip.background = element_blank(),
      strip.text       = element_text(size = 12),
      axis.text.x      = element_text(size = 10, angle = 30, hjust = 1),
      axis.text.y      = element_text(size = 10),
      axis.title       = element_text(size = 14),
      legend.position  = "none",
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    )
}

# Arabidopsis genes annotated with GO:0009853 photorespiration
photo_canonical <- c("AT5G36700","AT5G47760",                          # PGLP
                     "AT3G14415","AT3G14420","AT4G18360",              # GOX
                     "AT1G23310","AT1G70580","AT2G13360",              # GGT
                     "AT4G33010","AT2G26080","AT1G11860","AT2G35370",
                     "AT1G32470","AT3G17240",                          # GDC
                     "AT4G37930","AT5G26780",                          # SHMT
                     "AT1G68010","AT1G79870",                          # HPR
                     "AT1G80380",                                      # GLYK
                     "AT1G67090","AT5G38410","AT5G38420","AT5G38430")  # RBCS

photo_family <- c("AT5G36700"="PGLP","AT5G47760"="PGLP",
                  "AT3G14415"="GOX","AT3G14420"="GOX","AT4G18360"="GOX",
                  "AT1G23310"="GGT","AT1G70580"="GGT","AT2G13360"="GGT",
                  "AT4G33010"="GDC","AT2G26080"="GDC","AT1G11860"="GDC",
                  "AT2G35370"="GDC","AT1G32470"="GDC","AT3G17240"="GDC",
                  "AT4G37930"="SHMT","AT5G26780"="SHMT",
                  "AT1G68010"="HPR","AT1G79870"="HPR",
                  "AT1G80380"="GLYK",
                  "AT1G67090"="RBCS","AT5G38410"="RBCS",
                  "AT5G38420"="RBCS","AT5G38430"="RBCS")

# Arabidopsis genes annotated with GO:0010011 auxin binding
auxin_canonical <- c("AT3G62980","AT4G03190","AT3G26810","AT1G12820",
                     "AT4G24390","AT5G49980",                          # TIR/AFB
                     "AT2G38120","AT2G21050","AT1G77690","AT1G54990",  # AUX/LAX
                     "AT4G02980")                                      # ABP1

auxin_family <- c("AT3G62980"="TIR/AFB","AT4G03190"="TIR/AFB",
                  "AT3G26810"="TIR/AFB","AT1G12820"="TIR/AFB",
                  "AT4G24390"="TIR/AFB","AT5G49980"="TIR/AFB",
                  "AT2G38120"="AUX/LAX","AT2G21050"="AUX/LAX",
                  "AT1G77690"="AUX/LAX","AT1G54990"="AUX/LAX",
                  "AT4G02980"="ABP1")

p_photo_fam <- core_family_plot(
  photo_canonical, photo_family,
  c("PGLP","GOX","GGT","GDC","SHMT","HPR","GLYK","RBCS"),
  "photorespiration", "Q2a_up_in_hot",
  c("FALSE" = "#C56B8A", "TRUE" = "#2C5F8D"),
  "Core Photorespiration Homologs", 4)

p_auxin_fam <- core_family_plot(
  auxin_canonical, auxin_family,
  c("TIR/AFB","AUX/LAX","ABP1"),
  "auxin binding", "Q2a_up_in_cold",
  c("FALSE" = "#5A9367", "TRUE" = "#2C5F8D"),
  "Core Auxin Binding Homologs", 3)

pdf("figS13a_photorespiration_core.pdf", width = 6, height = 4.5,
    useDingbats = FALSE, bg = "transparent")
print(p_photo_fam); dev.off()

pdf("figS13b_auxin_core.pdf", width = 4.5, height = 3.5,
    useDingbats = FALSE, bg = "transparent")
print(p_auxin_fam); dev.off()