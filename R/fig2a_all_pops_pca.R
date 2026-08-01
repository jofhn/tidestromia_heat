# Figure 2A-B — species-range PCA by collection region and by survival group
# Inputs: data/all_pops_imputed_speciesrange_ne_150kb.eigenvec
#         data/all_pops_imputed_speciesrange_ne_150kb.eigenval
#         data/sample_metadata.csv   (id, region, group)
#
# Eigenvectors are from PLINK 2 (see wgs/survival_panel.md, section 11).

library(data.table)
library(ggplot2)

PCX <- 1
PCY <- 2

region_cols <- c(DV = "#EE6677", perimeter = "#4477AA", outside = "#228833")
group_cols  <- c(survivor = "#457B9D", early_nonsurvivor = "#E63946",
                 late_nonsurvivor = "#F4A261")

# --- data --------------------------------------------------------------------
eigvec <- fread("data/all_pops_imputed_speciesrange_ne_150kb.eigenvec")
setnames(eigvec, "#IID", "id")
eigval <- scan("data/all_pops_imputed_speciesrange_ne_150kb.eigenval", quiet = TRUE)

# library names are <run>_<id>_S<n>; keep the numeric sample id
eigvec$id <- as.character(as.integer(sub("^.*_(\\d+)_S\\d+$", "\\1", eigvec$id)))

samples <- read.csv("data/sample_metadata.csv", colClasses = "character")

pca_df <- data.frame(
  PCX    = eigvec[[paste0("PC", PCX)]],
  PCY    = eigvec[[paste0("PC", PCY)]],
  region = samples$region[match(eigvec$id, samples$id)],
  group  = samples$group[match(eigvec$id, samples$id)]
)
pca_df <- pca_df[!is.na(pca_df$region), ]
pca_df$region <- factor(pca_df$region, levels = names(region_cols))
pca_df$group  <- factor(pca_df$group,
                        levels = c("survivor", "spacer",
                                   "early_nonsurvivor", "late_nonsurvivor"))

# shuffle so no category is drawn consistently on top
set.seed(42)
pca_df <- pca_df[sample(nrow(pca_df)), ]

xl <- paste0("PC", PCX, " (", round(eigval[PCX] / sum(eigval) * 100, 1), "%)")
yl <- paste0("PC", PCY, " (", round(eigval[PCY] / sum(eigval) * 100, 1), "%)")

pca_theme <- theme_classic(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position      = "top",
    legend.justification = "left",
    legend.text          = element_text(size = 12, margin = margin(l = -3)),
    legend.key           = element_blank(),
    legend.key.spacing.x = unit(-0.1, "pt"),
    legend.box.spacing   = unit(2, "pt"),
    legend.margin        = margin(l = -30),
    plot.margin          = margin(r = 10),
    axis.text            = element_text(size = 10),
    axis.title           = element_text(size = 14),
    panel.background     = element_rect(fill = "transparent", color = NA),
    plot.background      = element_rect(fill = "transparent", color = NA),
    legend.background    = element_rect(fill = "transparent", color = NA)
  )

# --- 2A: by collection region ------------------------------------------------
p_region <- ggplot(pca_df, aes(PCX, PCY, color = region)) +
  geom_point(size = 2.5, alpha = 0.8, shape = 16) +
  scale_color_manual(values = region_cols, name = "",
                     labels = c("DV", "Perimeter", "Outside"),
                     breaks = names(region_cols)) +
  labs(x = xl, y = yl) +
  pca_theme

pdf("fig2a_all_pops_pca_region.pdf", width = 2.75, height = 2.75, useDingbats = FALSE)
print(p_region)
dev.off()

# --- 2B: by survival group ---------------------------------------------------
# empty "spacer" level fills the unused slot in the two-column legend
p_group <- ggplot(pca_df, aes(PCX, PCY, color = group)) +
  geom_point(size = 2.5, alpha = 0.8, shape = 16) +
  scale_color_manual(
    values = c(group_cols, spacer = NA),
    name   = "",
    labels = c(survivor = "Survivor", spacer = "",
               early_nonsurvivor = "Early Nonsurvivor",
               late_nonsurvivor  = "Late Nonsurvivor"),
    breaks = c("survivor", "spacer", "early_nonsurvivor", "late_nonsurvivor"),
    drop   = FALSE
  ) +
  guides(color = guide_legend(ncol = 2, byrow = FALSE)) +
  labs(x = xl, y = yl) +
  pca_theme +
  theme(legend.key.spacing.y = unit(-5, "pt"))

pdf("fig2a_all_pops_pca_group.pdf", width = 2.75, height = 2.75, useDingbats = FALSE)
print(p_group)
dev.off()
