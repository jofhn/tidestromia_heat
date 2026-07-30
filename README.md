# tidestromia_heat

Code for the analyses in Feehan et al. (bioRxiv 2026.07.22.740203), on heat
adaptation in *Tidestromia oblongifolia*.

## WGS

Raw WGS reads will be released at NCBI SRA upon publication.
The AmA10 assembly and annotation are from Prado et al. 2025 (Curr Biol 35:5502-5520).

`wgs/reference_panel.md` covers the reference panel resequencing: read
processing, variant calling and filtering, then PCA and F_ST.

`wgs/survival_panel.md` covers the survival panel resequencing: read
processing, variant calling and filtering, imputation against the reference panel,
PCA and F_ST, and genome-wide association testing.

## RNA-seq

Raw reads and count matrices will be released at NCBI GEO upon publication.

`rnaseq/hotcold.md` covers the heat/cold series: read processing, STAR genome
index and alignment, and assembly of the gene-by-sample count matrix.

## R-scripts

### Required packages

Scripts were run under R 4.3.2. Figures use Helvetica.

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "data.table", "readxl", "stringr",
  "sf", "terra", "tidyterra", "geodata", "rnaturalearth", "rnaturalearthdata",
  "ggnewscale", "ggrepel", "ggtext", "ggh4x", "glue", "eulerr", "circlize",
  "rgbif", "rinat", "httr", "jsonlite",
  "lmerTest", "emmeans", "multcompView", "FSA", "MCMCglmm",
  "survival", "survminer", "ashr", "ragg", "magick"
))

# not on CRAN; needed by ne_states()
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")

# Bioconductor
BiocManager::install(c("DESeq2", "apeglm", "RUVSeq", "clusterProfiler",
                       "ComplexHeatmap"))
```

`sf` and `terra` need system GDAL, GEOS and PROJ. `eulerr` >= 8.0 requires Rust
(https://rustup.rs).

### Figure scripts

`R/` holds one script per main-text figure panel. Input CSVs are in `R/data/`;
scripts that fetch their data at runtime need no input file.

`R/fig1a_range_map.R` — Fig. 1A. GBIF and iNaturalist occurrence records over
WorldClim summer maximum temperature. All inputs are fetched at runtime.
WorldClim rasters are cached in the working directory on first run.

`R/fig1b_collection_map.R` — Fig. 1B. Collection sites over shaded terrain,
coloured by collection region. Reads `data/collection_regions.csv`.

`R/fig1c_collection_curves.R` — Fig. 1C. Mean hourly summer temperature by
collection region for June, July and August, with a mixed model over daytime
hours. Reads the output of `R/fetch_era_land5.R`.

`R/fig1d_collection_boxplot.R` — Fig. 1D. Daily maximum temperature by
collection region and month, with compact letter display from estimated
marginal means. Reads the output of `R/fetch_era_land5.R`.

`R/fig1g_range_map.R` — Fig. 1G. Health status scores across the lethal heat
trial, normalised within each day. Reads
`data/fig1g_survival_health_status_scores.csv`.

`R/fig2a_all_pops_pca.R` — Fig. 2A. Species-range PCA coloured by collection
region (left) and survival group (right). Reads the PLINK eigenvectors and
`data/sample_metadata.csv`.

`R/fig2b_qq.R` — Fig. 2B. QQ plot of the DV/Perimeter association test, with
genomic inflation factor and the Bonferroni threshold over LD-pruned variants.
The association results are not deposited (1.4 GB); see `wgs/survival_panel.md`.

`R/fig2c_manhattan.R` — Fig. 2C. Manhattan plot of the DV/Perimeter association
test, with clumped regions and lead SNPs highlighted.

`R/fig2d_loci.R` — Fig. 2D. LD decay around each lead SNP with candidate gene
positions, and the LD block boundaries derived from it.

`R/fig3a_survival.R` — Fig. 3A. Survival curves by thermal cluster. Reads
`data/tray1_tray2_tidy_temps.csv`.

`R/fig3bc_IR_clusters.R` — Fig. 3B-C. Leaf temperature depression by thermal
cluster at 60 °C, and the pairwise cluster contrasts.

`R/fig3d_MCMCglmm.R` — Fig. 3D. Bayesian generalised mixed model of the
predictors of death at 60 °C.

`R/fig3e_substrate.R` — Fig. 3E. Leaf temperature depression of survivors
against bare substrate, as a time series and as per-plant extremes.

`R/fig4_RNAseq.R` — Fig. 4 and figs. S12-S13. Differential expression across
the heat series: sample PCA, the 55 °C heatmap, GO enrichment, and expression
trajectories of enriched and core pathway gene sets.

`R/fig5a_beta.R` — Fig. 5A. Leaf temperature depression across air temperatures
by thermal cluster, and the per-plant slope of that depression against air
temperature.

`R/fig5b_energy_balance.R` — Fig. 5B. Leaf energy balance model: leaf-air
temperature difference and transpiration across stomatal conductance. No input
file; the model is parameterised in the script.

`R/fig5c_ABA.R` — Fig. 5C. Leaf temperature depression under ABA, mock and
untreated painting at 60 °C, as a time series and as per-plant extremes.

`R/fig5d_ABA_survival.R` — Fig. 5D. Difference in survival between the painted
leaf and its internal control, by treatment.