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

`R/` holds one script per main-text figure panel. Input CSVs are in `R/data/`;
scripts that fetch their data at runtime need no input file.

`R/fig1a_range_map.R` — Fig. 1A. GBIF and iNaturalist occurrence records over
WorldClim summer maximum temperature. Fetches all inputs; no data file.


