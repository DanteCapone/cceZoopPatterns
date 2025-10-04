# cceZoopPatterns
Zooplankton community patterns in the California Current Ecosystem from eDNA metabarcoding across environmental gradients and size fractions.

## Summary
This repository contains the analysis code and minimal input data for a manuscript on spatial and environmental patterns in zooplankton communities of the California Current Ecosystem (CCE) using eDNA metabarcoding. We integrate 18S rRNA and COI markers, environmental metadata (e.g., distance from shore, nitracline depth, hypoxia depth, chl max depth), and phytoplankton covariates. Analyses are performed by size fraction (S1: 0.2–0.5 mm, S2: 0.5–1 mm, S3: 1–2 mm), leveraging compositional modeling (fido/pibble with CLR) and standard community ecology workflows (phyloseq, vegan). We compare metabarcoding-derived patterns with Zooscan-derived biomass and visualize key taxa–environment relationships (heatmaps, effect sizes) and community composition (treemaps).

## What’s included (no FASTA/FASTQ)
We include only the minimal CSV inputs needed to run the main analyses. Raw sequencing files (FASTA/FASTQ) and private data are excluded.

- inst/extdata/physical_environmental_data/
  - env_metadata_impute_phyloseq_6.9.2023.csv
  - env_metadata_impute_phyloseq_6.2.2023_for_map.csv
  - sample_depths.csv
  - monthly_composite_predicted_PC1.csv
- inst/extdata/raw_reads/
  - ASV_table_18s_run1.csv(.gz), ASV_table_18s_run2.csv(.gz)
  - ASV_table_coi_run1.csv(.gz), ASV_table_coi_run2.csv(.gz)
- inst/extdata/taxa_files/
  - blast_metazoo_18s_holo.csv
  - blast_metazoo_coi_holo.csv
- inst/extdata/fido/phy/
  - fido_18s_s{1,2,3}_family_phy_all_subpools(_nofilt).csv (as used by scripts)
  - fido_coi_s{1,2,3}_ecdf_taxa_phy.csv (as used by scripts)
- inst/extdata/Zooscan/
  - zooscan_by_sample_biomass_esd.csv

Note: If any file approaches 100 MB, we store a compressed `.csv.gz`. The analysis code transparently reads gzipped CSVs.

## How to reproduce
1) Clone and open the project root.
2) Start R in the project root (R >= 4.2 recommended).
3) Restore dependencies:
```r
install.packages("renv")
renv::init()          # first time
renv::restore()       # thereafter
install.packages("devtools")
devtools::load_all()  # expose functions in R/
```
4) Run the main analysis:
```r
source("scripts/Zoop_Patterns_main_analysis.R")
```
Outputs are written to:
- figures/
- data/processed/

## Repository structure
- R/ … exported helper functions (formerly in scripts/helpful_functions/)
- scripts/
  - Zoop_Patterns_main_analysis.R (entry point)
- inst/extdata/ … small, versioned inputs used by scripts
- data/
  - raw/ (gitignored)
  - processed/ (generated outputs)
- figures/ (generated outputs)
- renv/, renv.lock … locked dependency versions
- DESCRIPTION, NAMESPACE, .Rbuildignore, .Rprofile, .gitignore

## Methods (brief)
- Data integration and wrangling via tidyverse.
- Community objects and transforms via phyloseq and microbiome (CLR).
- Taxa–environment modeling via fido (pibble; Gamma hyperparameter scanning; posterior summaries).
- Visualization via ggplot2, ggdist, gridExtra; community composition via treemaps.
- Size-fraction analyses (S1–S3) and comparisons to Zooscan biomass (ESD-binned).

## Data policy
- Raw FASTA/FASTQ, BAM/SAM, and private data are excluded.
- Only small derivative inputs are versioned under inst/extdata/.
- Scripts write to data/processed/ and must never overwrite data/raw/ (see .gitignore).

## Citation
- Paper: “Zooplankton community patterns in the CCE from eDNA metabarcoding” (in prep).
- Please cite this repository and the paper once published. A CITATION.cff will be added upon acceptance.

## License
MIT License (see LICENSE).

## Contact
Questions: open a GitHub issue or contact the corresponding author(s) listed in the manuscript.
