# Fido Analysis: NBSS Slope and Top SIMPER Taxa
# Uses phyloseq objects, top taxa from capscale SIMPER, and NBSS slope as covariate
# Follows original model structure: split by size, merge by Sample_ID_short
# Author: Dante Capone
# Date: 2024

# Load packages -----------------------------------------------------------
librarian::shelf(
  tidyverse,
  here,
  phyloseq,
  fido,
  microbiome,
  patchwork
)

base_dir <- here("scripts/NBSS_fido_analysis")
data_dir <- file.path(base_dir, "data")
outputs_dir <- file.path(base_dir, "outputs")
figures_dir <- file.path(base_dir, "figures")

# Load helper functions
source(file.path(base_dir, "fido_missing_funs.R"))

# Load data ---------------------------------------------------------------

# Load phyloseq objects
Phy_coi_raw <- readRDS(file.path(data_dir, "ps_raw_coi.rds"))
Phy_18s_raw <- readRDS(file.path(data_dir, "ps_raw_18s.rds"))

cat("Loaded COI phyloseq:", ntaxa(Phy_coi_raw), "taxa,", nsamples(Phy_coi_raw), "samples\n")
cat("Loaded 18S phyloseq:", ntaxa(Phy_18s_raw), "taxa,", nsamples(Phy_18s_raw), "samples\n")

# Load NBSS slopes
nbss_slopes <- read.csv(file.path(data_dir, "normalized_biomass_spectrum_slopes_esd.csv"))
if (interactive()) view(nbss_slopes)
cat("Loaded NBSS slopes for", nrow(nbss_slopes), "samples\n")

# Load top SIMPER taxa from capscale analysis
simper_top_taxa_coi <- read.csv(file.path(data_dir, "simper_top_taxa_clusters.csv"))
simper_top_taxa_18s <- read.csv(file.path(data_dir, "simper_top_taxa_clusters_18S.csv"))

cat("Loaded", nrow(simper_top_taxa_coi), "COI top SIMPER taxa\n")
cat("Loaded", nrow(simper_top_taxa_18s), "18S top SIMPER taxa\n")

# Prepare NBSS slopes for merging -----------------------------------------

# Standardize sample IDs for joining with NBSS (NBSS is at Sample_ID_short level)
nbss_for_merge <- nbss_slopes %>%
  select(Sample_ID_short, slope)

cat("Prepared NBSS slopes for", nrow(nbss_for_merge), "unique Sample_ID_short values\n")

# Extract and prepare top SIMPER taxa for COI (species level) -------------

top_species_coi <- simper_top_taxa_coi %>%
  pull(species_label) %>%
  unique() %>%
  na.omit() %>%
  as.character()

cat("Found", length(top_species_coi), "unique COI species from SIMPER\n")

# Extract and prepare top SIMPER taxa for 18S (genus level) ---------------

top_genera_18s <- simper_top_taxa_18s %>%
  pull(genus_label) %>%
  unique() %>%
  na.omit() %>%
  as.character()

cat("Found", length(top_genera_18s), "unique 18S genera from SIMPER\n")

# Function to fix taxa names and subset to SIMPER taxa --------------------

prepare_phyloseq_for_fido <- function(ps_obj, top_taxa_list, marker_name = "COI", 
                                      glom_rank = "Species") {
  
  cat("\n=== Preparing", marker_name, "phyloseq ===\n")
  
  # Extract taxonomy table
  tax_df <- as.data.frame(tax_table(ps_obj)) %>%
    rownames_to_column("ASV")
  
  # Prepare taxa for matching based on glom_rank
  if (glom_rank == "Species") {
    # For COI: Fix specific SIMPER taxa that need renaming (Acartia and Sagitta)
    tax_df <- tax_df %>%
      mutate(
        Taxa_fixed = case_when(
          # If Species is NA or contains "Unidentified" and Genus is Acartia or Sagitta, rename
          (is.na(Species) | grepl("Unidentified", Species, ignore.case = TRUE)) & Genus == "Acartia" ~ "Unidentified Acartia",
          (is.na(Species) | grepl("Unidentified", Species, ignore.case = TRUE)) & Genus == "Sagitta" ~ "Unidentified Sagitta",
          # Otherwise use Species as is
          TRUE ~ Species
        )
      )
  } else if (glom_rank == "Genus") {
    # For 18S: Use Genus directly
    tax_df <- tax_df %>%
      mutate(Taxa_fixed = Genus)
  }
  
  # Filter to taxa matching SIMPER list
  matching_asvs <- tax_df %>%
    filter(Taxa_fixed %in% top_taxa_list) %>%
    pull(ASV)
  
  cat("Found", length(matching_asvs), "taxa matching SIMPER", glom_rank, "\n")
  
  if (length(matching_asvs) == 0) {
    warning(paste("No matching taxa found for", marker_name))
    return(NULL)
  }
  
  # Subset phyloseq to matching taxa
  ps_subset <- prune_taxa(matching_asvs, ps_obj)
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
  
  cat("Before glomming:", ntaxa(ps_subset), "taxa\n")
  
  # Glom at specified rank
  ps_subset <- tax_glom(ps_subset, taxrank = glom_rank)
  
  cat("After", glom_rank, "-level glomming:", ntaxa(ps_subset), "taxa\n")
  
  # Get the taxa names for the glommed taxa and make them unique
  tax_table_subset <- as.data.frame(tax_table(ps_subset))
  
  # Create final names based on glom_rank
  if (glom_rank == "Species") {
    taxa_names_final <- tax_table_subset %>%
      mutate(
        Taxa_fixed = case_when(
          (is.na(Species) | grepl("Unidentified", Species, ignore.case = TRUE)) & Genus == "Acartia" ~ "Unidentified Acartia",
          (is.na(Species) | grepl("Unidentified", Species, ignore.case = TRUE)) & Genus == "Sagitta" ~ "Unidentified Sagitta",
          TRUE ~ Species
        )
      ) %>%
      pull(Taxa_fixed) %>%
      make.unique()
  } else {
    taxa_names_final <- tax_table_subset %>%
      pull(Genus) %>%
      make.unique()
  }
  
  # Replace taxa_names (rownames) with fixed names
  taxa_names(ps_subset) <- taxa_names_final
  
  cat("Final subset:", ntaxa(ps_subset), "taxa with", glom_rank, "-level names\n")
  cat("Taxa names:", paste(taxa_names(ps_subset), collapse = ", "), "\n")
  
  # Return phyloseq with named vector of taxa names
  fixed_names <- setNames(taxa_names_final, taxa_names_final)
  
  return(list(ps = ps_subset, taxa_names = fixed_names))
}

# Prepare both markers
coi_prep <- prepare_phyloseq_for_fido(Phy_coi_raw, top_species_coi, 
                                      marker_name = "COI", glom_rank = "Species")
s18_prep <- prepare_phyloseq_for_fido(Phy_18s_raw, top_genera_18s, 
                                      marker_name = "18S", glom_rank = "Genus")

Phy_coi_subset <- coi_prep$ps
Phy_18s_subset <- s18_prep$ps

# Split by size fraction and prepare for fido ----------------------------

# Function to split by size and run fido models
run_fido_by_size <- function(ps_obj, marker_name = "COI") {
  
  cat("\n=== Running Fido models for", marker_name, "by size fraction ===\n")
  
  # Split by size fraction (S1, S2, S3)
  data_list <- list(
    S1 = prune_samples(grepl("S1", sample_names(ps_obj)), ps_obj),
    S2 = prune_samples(grepl("S2", sample_names(ps_obj)), ps_obj),
    S3 = prune_samples(grepl("S3", sample_names(ps_obj)), ps_obj)
  )
  
  # Merge replicates by Sample_ID_short using mean
  data_list <- lapply(data_list, function(ps) {
    
    # Get unique Sample_ID_short values before merging
    unique_samples <- unique(sample_data(ps)$Sample_ID_short)
    
    # Merge samples
    ps_merged <- merge_samples(ps, "Sample_ID_short", fun = mean)
    
    # Re-add metadata with NBSS slopes
    meta_df <- data.frame(
      Sample_ID_short = sample_names(ps_merged),
      row.names = sample_names(ps_merged)
    ) %>%
      left_join(nbss_for_merge, by = "Sample_ID_short") %>%
      column_to_rownames("Sample_ID_short")
    
    sample_data(ps_merged) <- sample_data(meta_df)
    
    return(ps_merged)
  })
  
  set.seed(899)
  
  # Loop through fido models for each size
  posterior_summaries <- list()
  
  for (i in seq_along(data_list)) {
    
    dat_name <- names(data_list)[i]
    cat("\n--- Size fraction:", dat_name, "---\n")
    
    dat <- data_list[[dat_name]]
    
    # Check if we have taxa
    if (ntaxa(dat) <= 1) {
      warning(paste("Not enough taxa for", marker_name, dat_name))
      next
    }
    
    # Prepare sample data with NBSS slope
    sample_dat <- as.data.frame(as(sample_data(dat), "matrix")) %>%
      filter(!is.na(slope))  # Only samples with NBSS data
    
    if (nrow(sample_dat) < 3) {
      warning(paste("Not enough samples with NBSS data for", marker_name, dat_name))
      next
    }
    
    cat("Using", nrow(sample_dat), "samples with NBSS slope data\n")
    
    # Select only slope as covariate
    sample_dat <- sample_dat %>%
      select(slope)
    
    # Build design matrix
    formula_obj <- as.formula(~ slope)
    X <- t(model.matrix(formula_obj, data = sample_dat))
    
    # Get OTU table
    Y <- otu_table(dat)
    if (!taxa_are_rows(dat)) {
      Y <- t(Y)
    }
    Y <- as.matrix(Y)
    
    # Match samples
    common_samples <- intersect(colnames(X), colnames(Y))
    X <- X[, common_samples, drop = FALSE]
    Y <- Y[, common_samples, drop = FALSE]
    
    cat("Design matrix:", nrow(X), "predictors x", ncol(X), "samples\n")
    cat("Count matrix:", nrow(Y), "taxa x", ncol(Y), "samples\n")
    
    # Use fixed gamma = 20 (as in original model)
    gamma <- 20
    
    # Prior specification
    upsilon <- ntaxa(dat) + 3
    Omega <- diag(ntaxa(dat))
    G <- cbind(diag(ntaxa(dat) - 1), -1)
    Xi <- (upsilon - ntaxa(dat)) * G %*% Omega %*% t(G)
    Theta <- matrix(0, ntaxa(dat) - 1, nrow(X))
    Gamma <- gamma * diag(nrow(X))
    
    # Fit model
    cat("Fitting pibble model...\n")
    priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
    priors <- to_clr(priors)
    
    names_covariates(priors) <- rownames(X)
    priors$Y <- Y
    
    posterior <- refit(priors, optim_method = "lbfgs", jitter = 1e-5)

    # Use taxa names from rownames of Y (already set as Species names)
    names_categories(posterior) <- rownames(Y)
    
    # Summarize posterior
    posterior_summary <- summary(posterior, pars = "Lambda")$Lambda
    
    # Add size column
    sizes <- c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")
    
    posterior_summary <- posterior_summary %>%
      mutate(size = sizes[i])
    
    posterior_summaries[[dat_name]] <- posterior_summary
    
    cat("Posterior summary complete for", dat_name, "\n")
  }
  
  # Combine all size fractions
  if (length(posterior_summaries) > 0) {
    posterior_summary_all <- bind_rows(posterior_summaries)
    
    # Extract taxonomy table for joining in plots
    tax_table_df <- as.data.frame(tax_table(ps_obj)) %>%
      rownames_to_column("taxa_name") %>%
      select(taxa_name, Order)
    
    return(list(
      posterior_summary_all = posterior_summary_all,
      posterior_summaries = posterior_summaries,
      taxonomy = tax_table_df
    ))
  } else {
    warning(paste("No successful models for", marker_name))
    return(NULL)
  }
}

# Run models for both markers ---------------------------------------------

# COI
results_coi <- run_fido_by_size(Phy_coi_subset, marker_name = "COI")

# 18S  
results_18s <- run_fido_by_size(Phy_18s_subset, marker_name = "18S")

# Visualize results -------------------------------------------------------

# Function to create coefficient plot with size fractions
plot_fido_results <- function(results, marker_name = "COI") {
  
  if (is.null(results)) {
    cat("No results to plot for", marker_name, "\n")
    return(NULL)
  }
  
  posterior_summary_all <- results$posterior_summary_all
  
  # Filter to significant results (95% CI doesn't include zero)
  sig_results <- posterior_summary_all %>%
    filter(sign(p2.5) == sign(p97.5)) %>%
    filter(covariate != "(Intercept)")
  
  if (nrow(sig_results) == 0) {
    cat("No significant results for", marker_name, "\n")
    return(NULL)
  }
  
  cat("Found", nrow(sig_results), "significant taxa-covariate associations for", marker_name, "\n")
  
  # Plot: grouped by taxon, colored by size
  p <- ggplot(sig_results, aes(x = mean, y = coord, color = size)) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.5, linewidth = 1.5) +
    geom_errorbarh(aes(xmin = p25, xmax = p75), height = 0.2, alpha = 0.8) +
    geom_point(aes(shape = size), size = 3, alpha = 0.9) +
    facet_wrap(~covariate, scales = "free_y", ncol = 1) +
    labs(
      title = paste(marker_name, "- NBSS Slope Effect on Top SIMPER Taxa"),
      x = "Centered Log-Ratio Change",
      y = "Taxon",
      color = "Size Fraction",
      shape = "Size Fraction"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.text.x = element_text(size = 9),
      strip.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
  
  return(p)
}

# Alternative: Plot by size fraction (faceted by size, all taxa on y-axis)
plot_fido_by_taxon <- function(results, marker_name = "COI") {
  
  if (is.null(results)) {
    return(NULL)
  }
  
  posterior_summary_all <- results$posterior_summary_all
  taxonomy <- results$taxonomy
  
  # Filter to slope covariate only (exclude intercept)
  sig_results <- posterior_summary_all %>%
    filter(covariate == "slope") %>%
    mutate(taxa_name = str_remove(coord, "^clr_")) %>%
    # Abbreviate long species names for better plotting
    mutate(taxa_name = case_when(
      taxa_name == "Pleuromamma abdominalis edentata" ~ "Pleuromamma a. edentata",
      TRUE ~ taxa_name
    ))
  
  if (nrow(sig_results) == 0) {
    return(NULL)
  }
  
  # Join with taxonomy to get Order
  sig_results <- sig_results %>%
    left_join(taxonomy, by = "taxa_name")
  
  # Sort taxa by mean value in largest size class (1-2 mm)
  taxa_order <- sig_results %>%
    filter(size == "1-2 mm") %>%
    arrange(desc(mean)) %>%
    pull(taxa_name)
  
  sig_results <- sig_results %>%
    mutate(taxa_name = factor(taxa_name, levels = taxa_order))
  
  # Create color palette specific to marker - map related taxa to same colors across markers
  if (marker_name == "COI") {
    # COI species color mapping
    taxa_color_map <- c(
      "Calanus pacificus" = "#1f77b4",        # Calanoida - matches Calanus
      "Calanoida" = "#aec7e8",               # General Calanoida
      "Neocalanus gracilis" = "#ff7f0e",     # Calanoida - matches Neocalanus
      "Pleuromamma a. edentata" = "#2ca02c", # Calanoida - matches Pleuromamma
      "Clausocalanus furcatus" = "#d62728", # Calanoida - matches Clausocalanus
      "Ctenocalanus vanus" = "#9467bd",      # Calanoida - matches Ctenocalanus
      "Candacia bipinnata" = "#8c564b",      # Calanoida - matches Candacia
      "Nanomia bijuga" = "#e377c2",          # Siphonophorae
      "Eucalanus californicus" = "#7f7f7f",  # Calanoida - matches Eucalanus
      "Lensia campanella" = "#bcbd22",       # Siphonophorae
      "Nematoscelis difficilis" = "#17becf", # Euphausiacea
      "Euchaeta marina" = "#ff9896",         # Calanoida
      "Euphausia pacifica" = "#98df8a",      # Euphausiacea
      "Rosacea cymbiformis" = "#c5b0d5",     # Siphonophorae
      "Paraeuchaeta rubra" = "#c49c94",      # Calanoida
      "Mecynocera clausi" = "#f7b6d2",        # Calanoida
      "Ditrichocorycaeus anglicus" = "#dbdb8d", # Cyclopoida
      "Acrocalanus monachus" = "#9edae5"     # Calanoida
    )
    taxa_colors <- taxa_color_map[sig_results$taxa_name]
    names(taxa_colors) <- sig_results$taxa_name
  } else {
    # 18S genera color mapping  
    taxa_color_map <- c(
      "Mnemiopsis" = "#ffbb78",              # Lobata
      "Ctenocalanus" = "#9467bd",            # Calanoida - matches Ctenocalanus vanus
      "Rhincalanus" = "#ff9896",             # Calanoida
      "Clausocalanus" = "#d62728",           # Calanoida - matches Clausocalanus furcatus
      "Eucalanus" = "#7f7f7f",               # Calanoida - matches Eucalanus californicus
      "Calanus" = "#1f77b4",                 # Calanoida - matches Calanus pacificus
      "Calocalanus" = "#aec7e8",             # Calanoida
      "Paracalanus" = "#c49c94",             # Calanoida
      "Metridia" = "#f7b6d2",                # Calanoida
      "Pleuromamma" = "#2ca02c",             # Calanoida - matches Pleuromamma abdominalis edentata
      "Scolecithricella" = "#9edae5",        # Calanoida
      "Neocalanus" = "#ff7f0e",              # Calanoida - matches Neocalanus gracilis
      "Candacia" = "#8c564b",                # Calanoida - matches Candacia bipinnata
      "Oithona" = "#dbdb8d",                 # Cyclopoida
      "Tomopteris" = "#17becf"               # Phyllodocida
    )
    taxa_colors <- taxa_color_map[sig_results$taxa_name]
    names(taxa_colors) <- sig_results$taxa_name
  }
  
  # Plot: x = taxa (sorted), y = mean, shape by size fraction, colored by taxon
  p <- ggplot(sig_results, aes(y = mean, x = taxa_name, color = taxa_name, fill = taxa_name, shape = size)) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5, linewidth = 1) +
    # 95% credible interval (p2.5 to p97.5) - lighter shading
    geom_linerange(aes(ymin = p2.5, ymax = p97.5, group = size, color = taxa_name),
                   linewidth = 6, alpha = 0.2,
                   position = position_dodge(width = 0.8)) +
    # 50% credible interval (p25 to p75) - darker shading, same color as point
    geom_linerange(aes(ymin = p25, ymax = p75, group = size, color = taxa_name),
                   linewidth = 6, alpha = 0.5,
                   position = position_dodge(width = 0.8)) +
    geom_point(size = 4, alpha = 1, color = "black", stroke = 0.8, position = position_dodge(width = 0.8)) +
    scale_color_manual(values = taxa_colors, guide = "none") +
    scale_fill_manual(values = taxa_colors, guide = "none") +
    scale_shape_manual(values = c(21, 24, 22), name = "Size Fraction") +
    labs(
      title = paste(marker_name),
      y = "Centered Log-Ratio with NBSS slope",
      x = "Taxon"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(size = 8, angle = 30, hjust = 1, face = "italic"),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold")
    )
  
  return(p)
}

# Create plots for both markers
if (!is.null(results_coi)) {
  plot_coi <- plot_fido_results(results_coi, "COI")
  if (!is.null(plot_coi)) print(plot_coi)
  
  plot_coi_taxon <- plot_fido_by_taxon(results_coi, "COI")
  if (!is.null(plot_coi_taxon)) print(plot_coi_taxon)
}

if (!is.null(results_18s)) {
  plot_18s <- plot_fido_results(results_18s, "18S")
  if (!is.null(plot_18s)) print(plot_18s)
  
  plot_18s_taxon <- plot_fido_by_taxon(results_18s, "18S")
  if (!is.null(plot_18s_taxon)) print(plot_18s_taxon)
}

# Combine plots with patchwork --------------------------------------------------

# Create combined plot with 2 rows (COI top, 18S bottom)
if (exists("plot_coi_taxon") && exists("plot_18s_taxon") && 
    !is.null(plot_coi_taxon) && !is.null(plot_18s_taxon)) {
  
  combined_fido_plot <- (plot_coi_taxon / plot_18s_taxon) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  print(combined_fido_plot)
  
  # Save combined plot
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_combined_taxa.pdf"),
    plot = combined_fido_plot,
    width = 12,
    height = 10,
    units = "in"
  )
  
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_combined_taxa.png"),
    plot = combined_fido_plot,
    width = 12,
    height = 10,
    units = "in",
    dpi = 600
  )
  
  cat("Saved combined COI/18S FIDO plot\n")
}

# Save results ------------------------------------------------------------

# Create output directories if they don't exist
dir.create(file.path(outputs_dir, "fido"), showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Save posterior summaries
if (!is.null(results_coi)) {
  write.csv(
    results_coi$posterior_summary_all,
    file.path(outputs_dir, "fido", "fido_nbss_coi_posterior_summary.csv"),
    row.names = FALSE
  )
  cat("Saved COI posterior summary\n")
}

if (!is.null(results_18s)) {
  write.csv(
    results_18s$posterior_summary_all,
    file.path(outputs_dir, "fido", "fido_nbss_18s_posterior_summary.csv"),
    row.names = FALSE
  )
  cat("Saved 18S posterior summary\n")
}

# Save individual plots
if (!is.null(results_coi) && !is.null(plot_coi)) {
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_coi.pdf"),
    plot = plot_coi,
    width = 8,
    height = 6,
    units = "in"
  )
  
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_coi.png"),
    plot = plot_coi,
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
  )
}

if (!is.null(results_coi) && !is.null(plot_coi_taxon)) {
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_coi_by_taxon.pdf"),
    plot = plot_coi_taxon,
    width = 10,
    height = 8,
    units = "in"
  )
  
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_coi_by_taxon.png"),
    plot = plot_coi_taxon,
    width = 10,
    height = 8,
    units = "in",
    dpi = 600
  )
}

if (!is.null(results_18s) && !is.null(plot_18s)) {
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_18s.pdf"),
    plot = plot_18s,
    width = 8,
    height = 6,
    units = "in"
  )
  
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_18s.png"),
    plot = plot_18s,
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
  )
}

if (!is.null(results_18s) && !is.null(plot_18s_taxon)) {
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_18s_by_taxon.pdf"),
    plot = plot_18s_taxon,
    width = 10,
    height = 8,
    units = "in"
  )
  
  ggsave(
    filename = file.path(figures_dir, "fido_nbss_18s_by_taxon.png"),
    plot = plot_18s_taxon,
    width = 10,
    height = 8,
    units = "in",
    dpi = 600
  )
}

cat("\n=== Fido analysis complete ===\n")
