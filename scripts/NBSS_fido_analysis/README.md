# NBSS + FIDO analysis (Zooscan + metabarcoding)

This folder contains a self-contained analysis that relates **Zooscan Normalized Biomass Size Spectrum (NBSS) slope** to **eDNA community composition**, using **FIDO** (multinomial logistic-normal / pibble models) on a curated set of **top SIMPER taxa**.

## What this analysis does

- **Input data**
  - Two phyloseq objects:
    - COI: `ps_raw_coi.rds`
    - 18S: `ps_raw_18s.rds`
  - NBSS slopes derived from Zooscan:
    - `normalized_biomass_spectrum_slopes_esd.csv` (uses the `slope` column)
  - Top taxa lists from the capscale/SIMPER workflow:
    - COI: `simper_top_taxa_clusters.csv` (uses `species_label`)
    - 18S: `simper_top_taxa_clusters_18S.csv` (uses `genus_label`)

- **Processing steps**
  - Extracts the taxa in the phyloseq objects that match the SIMPER top taxa list.
    - COI: works at **Species** level (with specific name-fixes for *Acartia* and *Sagitta* “Unidentified ...”).
    - 18S: works at **Genus** level.
  - Splits samples by size fraction (`S1`, `S2`, `S3`) using `sample_names()` pattern matching.
  - Merges technical replicates by `Sample_ID_short` (mean), then joins NBSS slope metadata.
  - Fits a FIDO pibble model with a single covariate: `~ slope`.

- **Outputs**
  - Posterior summaries (CSV):
    - `outputs/fido/fido_nbss_coi_posterior_summary.csv`
    - `outputs/fido/fido_nbss_18s_posterior_summary.csv`
  - Figures (PDF/PNG) saved to `figures/`.

## Folder structure

- `fido_nbss_top_taxa.R`
  - main analysis script
- `fido_missing_funs.R`
  - helper functions required for the FIDO workflow
- `data/`
  - required input files for running this script
- `outputs/`
  - model outputs written by the script
- `figures/`
  - figures written by the script

## How to run

1. Open / set your working directory to the repo root:
   - `cceZoopPatterns_repo`

2. Run the script:
   - `scripts/NBSS_fido_analysis/fido_nbss_top_taxa.R`

Notes:
- This script uses `here()`; it expects to be run with the repo root as the project root.
- FIDO model fitting can take time depending on the number of taxa retained and your machine.

## Provenance of the input files

The `data/` inputs in this folder are copies of the files used in the original CCE Zooplankton Metabarcoding publication repo, placed here to make this analysis easier to share/re-run.
