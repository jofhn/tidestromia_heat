# tidestromia_heat

Code for the analyses in Feehan et al. (bioRxiv 2026.07.22.740203), on heat
adaptation in *Tidestromia oblongifolia*.

Two input datasets too large for this repository are archived at Zenodo 
(https://doi.org/10.5281/zenodo.21729632):
The ERA5-Land temperature extraction (era5_land_temperature_2m.csv; Figs. 1B-D), and
the genome-wide association (GWA) testing results (glm_dvp_full.tsv.gz; Figs. 2B-C). 

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

`R/` holds scripts for figures. Input CSVs are in `R/data/`. 
Collection region and survival group assignments in `R/data/sample_metadata.csv`.

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

### Scripts

`R/fetch_era_land5.R` — fetches the ERA5-Land hourly temperature series for the
collection sites from the Open-Meteo archive API, or download the output from the
Zenodo deposit (https://doi.org/10.5281/zenodo.21729632).

`R/extract_leaf_temperature.R` — how leaf temperature was extracted from the
infrared snapshots, shown for one snapshot.

### Figure scripts

`R/fig1a_range_map.R` — Fig. 1A. GBIF and iNaturalist occurrence records over
WorldClim summer maximum temperature. All inputs are fetched at runtime.
WorldClim rasters are cached in the working directory on first run.

`R/fig1bcd_climate.R` — Figs. 1B-D. Collection sites over shaded terrain, mean
hourly summer temperature by collection region, and daily maximum temperature by
region and month. Requires the output of `R/fetch_era_land5.R` from Zenodo deposit
(https://doi.org/10.5281/zenodo.21729632).

`R/fig1g_survival_plot.R` — Fig. 1G. Extreme heat survival assay scoring. 

`R/fig2a_all_pops_pca.R` — Fig. 2A. Species-range PCA coloured by collection
region (left) and survival group (right).

`R/fig2b_qq.R` — Fig. 2B. QQ plot of the DV/Perimeter GWA model.
Requires `glm_dvp_full.tsv.gz` from Zenodo (https://doi.org/10.5281/zenodo.21729632)

`R/fig2c_manhattan.R` — Fig. 2C. Manhattan plot of the DV/Perimeter GWA model.
Requires glm_dvp_full.tsv.gz from Zenodo (https://doi.org/10.5281/zenodo.21729632)

`R/fig2d_loci.R` — Fig. 2D. LD decay around each lead SNP with candidate gene
positions, and the LD block boundaries derived from it.

`R/fig3a_survival_analysis.R` — Fig. 3A. Survival analysis by cluster.

`R/fig3bc_IR_clusters.R` — Fig. 3B-C. DeltaT by cluster at 60°C.

`R/fig3d_MCMCglmm.R` — Fig. 3D. Bayesian generalised mixed model of the
predictors of death at 60°C.

`R/fig3e_extended_cooling.R` — Fig. 3E. DeltaT of survivors after extended 60°C.

`R/fig4_RNAseq.R` — Fig. 4 and figs. S12-S13. Differential expression analysis.

`R/fig5a_beta.R` — Fig. 5A. DeltaT in changing air temperature.

`R/fig5b_energy_balance.R` — Fig. 5B. Leaf energy balance model.

`R/fig5c_ABA.R` — Fig. 5C. DeltaT during ABA treatment at 60°C.

`R/fig5d_ABA_survival.R` — Fig. 5D. Delta survival from ABA treatment.