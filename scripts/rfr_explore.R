# Load necessary libraries
library(tidyverse)
library(phyloseq)
library(randomForest)
library(caret)
library(here)

#Dryweights
biomass=read.csv(here("data/zoop_other/biomass_processed.csv")) %>% 
  select(c(Sample_ID_short, max_size, biomass_mg_m2))

# COI ---------------------------------------------------------------------



# Standardize OTU row names in leray_metazoo_otucoi and taxa_data
leray_metazoo_otucoi <- read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
rownames(leray_metazoo_otucoi) <- paste0("otu_", gsub(" ", ".", rownames(leray_metazoo_otucoi)))

# Standardize taxa table row names
leray_metazoo_taxa <- read.csv(here("data/phyloseq_bio_data/COI/coi_taxa_table_eDNA_metazoogene.csv")) %>%
  column_to_rownames("X")
rownames(leray_metazoo_taxa) <- paste0("otu_", gsub(" ", ".", rownames(leray_metazoo_taxa)))

# Standardize metadata row names
leray_metazoo_meta <- read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>% 
  left_join(.,biomass, by = c("Sample_ID_short","max_size"))%>%
  column_to_rownames("Sample_ID_dot") %>%
  dplyr::select(-X) )






# Convert to phyloseq object
OTU <- otu_table(as.matrix(leray_metazoo_otucoi), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(leray_metazoo_taxa))
meta <- sample_data(leray_metazoo_meta)
Phy_merged_coi <- phyloseq(OTU, TAX, meta)

# Extract OTU and taxa data
otu_data <- as(otu_table(Phy_merged_coi), "matrix")
taxa_data <- as.data.frame(tax_table(Phy_merged_coi))

# Replace spaces in taxa table columns
taxa_data <- taxa_data %>%
  mutate(across(everything(), ~ gsub(" ", ".", .)))  # Replace spaces with "."
meta_data <- as(sample_data(Phy_merged_coi), "data.frame") # Metadata

# # Step 3: Use species names as OTU table rownames
# species_names <- taxa_data$
# rownames(otu_data) <- species_names

# Normalize OTU table (relative abundance)
otu_normalized <- apply(otu_data, 2, function(x) x / sum(x))

# Step 4: Align with metadata and prepare for Random Forest
otu_samples <- colnames(otu_normalized)
meta_samples <- rownames(meta_data)
common_samples <- intersect(otu_samples, meta_samples)

otu_normalized <- otu_normalized[, common_samples]
meta_data <- meta_data[common_samples, ]

# Check for and handle missing values in PC1
if (any(is.na(meta_data$PC1))) {
  meta_data <- meta_data[!is.na(meta_data$PC1), ]
}

# Prepare data frame for modeling
rownames(otu_normalized) <- gsub(" ", ".", rownames(otu_normalized))
otu_df <- as.data.frame(t(otu_normalized)) # Transpose for predictors

# Ensure unique column names
colnames(otu_df) <- make.unique(colnames(otu_df))

# Merge columns containing "Species" using dplyr
otu_df <- otu_df %>%
  mutate(unidentified.spp = rowSums(select(., contains("Species"))))%>%
  select(-contains("Species"))

otu_df$PC1 <- meta_data$PC1

# Step 5: Train-Test Split
set.seed(123)
train_index <- createDataPartition(otu_df$PC1, p = 0.8, list = FALSE)
train_data <- otu_df[train_index, ]
test_data <- otu_df[-train_index, ]

# Step 6: Fit Random Forest Regression
set.seed(123)
rf_model <- randomForest(
  PC1 ~ ., 
  data = train_data,
  ntree = 500, 
  importance = TRUE
)

# Step 7: Evaluate Model
# Extract feature columns from training data
train_features <- train_data %>% select(-PC1)
test_features <- test_data %>% select(-PC1)

# Ensure columns match between training and test sets
if (!all(colnames(train_features) == colnames(test_features))) {
  stop("Column names in test data do not match the training data!")
}

# Predict on test data
predicted <- predict(rf_model, newdata = test_features)

# Compare predictions with actual values
actual <- test_data$PC1

# Calculate RMSE
rmse <- sqrt(mean((predicted - actual)^2))
cat("Root Mean Squared Error (RMSE):", rmse, "\n")

# Plot model performance
performance_plot <- data.frame(Actual = actual, Predicted = predicted)
ggplot(performance_plot, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Model Performance: Predicted vs Actual",
    x = "Actual PC1",
    y = "Predicted PC1"
  ) +
  theme_minimal()

# Step 8: Feature Importance
importance <- importance(rf_model)
importance_df <- data.frame(Feature = rownames(importance), Importance = importance[, 1]) %>%
  arrange(desc(Importance))

# Join importance_df with taxa_data to get species names
importance_with_taxa <- importance_df %>%
  left_join(taxa_data%>% rownames_to_column(var = "Feature"), by = c("Feature")) # Ensure taxa_data has row names or a matching column

# Replace spaces with "." in species names before handling duplicates
importance_with_taxa <- importance_with_taxa %>%
  mutate(Species = gsub(" ", ".", Species),
         Genus = gsub(" ", ".", Genus))%>%
# Replace Feature column with species names, or use genus names if species is missing
  mutate(Feature_label = ifelse(!is.na(Species), paste0(Species,"_",Feature), paste0("unidentified_",Genus,"_",Feature)))





# Filter out "unidentified_NA_NA" features
importance_with_taxa <- importance_with_taxa %>%
  filter(Feature != "unidentified_NA_NA")

# Order Feature by Importance for plotting
importance_with_taxa <- importance_with_taxa %>%
  mutate(Feature_label = factor(Feature_label, levels = Feature_label[order(Importance)]))

# Plot top 20 important features with species names
top_features <- head(importance_with_taxa, 20)
ggplot(top_features, aes(x = Feature_label, y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 Important Species Predicting PC1", x = "Species", y = "Importance") +
  theme_minimal()

# Optional: Extract the core species names or OTUs
core_otus <- top_features$Feature
cat("Core OTUs/Species:", paste(core_otus, collapse = ", "), "\n")



# Step 9: Loop through subsets of OTUs to find the best model
# Initialize storage for results
subset_results <- list()

# Rank OTUs by importance
ranked_otus <- importance_with_taxa$Feature

# Loop through decreasing numbers of OTUs
for (n_otus in seq(length(ranked_otus), 10, by = -10)) {
  
  # Subset top n OTUs
  selected_otus <- head(ranked_otus, n_otus)
  
  # Filter training and test datasets to include only the selected OTUs
  train_subset <- train_data %>% select(all_of(selected_otus), PC1)
  test_subset <- test_data %>% select(all_of(selected_otus), PC1)
  
  # Fit Random Forest model on the subset
  set.seed(123)
  rf_subset_model <- randomForest(
    PC1 ~ ., 
    data = train_subset,
    ntree = 500
  )
  
  # Predict and calculate RMSE
  subset_predicted <- predict(rf_subset_model, newdata = test_subset %>% select(-PC1))
  subset_actual <- test_subset$PC1
  subset_rmse <- sqrt(mean((subset_predicted - subset_actual)^2))
  
  
  # Store results
  subset_results[[as.character(n_otus)]] <- list(
    n_otus = n_otus,
    model = rf_subset_model,
    rmse = subset_rmse
  )
}

# Step 10: Summarize results and find the optimal model
results_df <- do.call(rbind, lapply(subset_results, function(res) {
  data.frame(n_otus = res$n_otus, RMSE = res$rmse)
}))

# Find the best model (min RMSE)
best_result <- results_df[which.min(results_df$RMSE), ]
cat("Best model uses", best_result$n_otus, "OTUs with RMSE:", best_result$RMSE, "\n")

# Plot RMSE vs. number of OTUs
ggplot(results_df, aes(x = n_otus, y = RMSE)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Model Performance vs Number of OTUs",
    x = "Number of OTUs",
    y = "RMSE"
  ) +
  theme_minimal()


# Evaluate the best model objectively
subset_rmse <- sqrt(mean((best_predicted - best_actual)^2))  # RMSE
subset_mae <- mean(abs(best_predicted - best_actual))        # MAE
subset_r2 <- 1 - (sum((best_actual - best_predicted)^2) / sum((best_actual - mean(best_actual))^2))  # R-squared

cat("Model Evaluation Metrics:\n")
cat("RMSE:", subset_rmse, "\n")
cat("MAE:", subset_mae, "\n")
cat("R-squared:", subset_r2, "\n")



# Find the best result (minimum RMSE)
best_result <- subset_results[[as.character(best_result$n_otus)]]

# Extract the best model and data
best_model <- best_result$model
best_test_subset <- test_data %>% select(all_of(head(ranked_otus, best_result$n_otus)), PC1)

# Predict using the best model
best_predicted <- predict(best_model, newdata = best_test_subset %>% select(-PC1))
best_actual <- best_test_subset$PC1

# Remake model performance plot
best_performance_plot <- data.frame(Actual = best_actual, Predicted = best_predicted)
ggplot(best_performance_plot, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(
    title = paste("Model Performance with Best Result (", best_result$n_otus, "OTUs)", sep = ""),
    x = "Actual PC1",
    y = "Predicted PC1"
  ) +
  theme_minimal()

# Get the top OTUs contributing to the best model
best_contributing_otus <- head(ranked_otus, best_result$n_otus)

# Filter importance_with_taxa for the top contributing OTUs
best_otu_importance <- importance_with_taxa %>%
  filter(Feature %in% best_contributing_otus)%>%
  mutate(Feature_label = substr(Feature_label, 1, 40))  # Truncate to first 20 characters

# Plot the importance of contributing OTUs
ggplot(best_otu_importance, aes(x = reorder(Feature_label, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = paste("Feature Importance for Best Model (", best_result$n_otus, "OTUs)", sep = ""),
    x = "OTUs/Species",
    y = "Importance"
  ) +
  theme_minimal()


# Define plot dimensions and resolution for L&O
plot_width <- 12  # Convert cm to inches (single-column width)
plot_height <- 4  # Adjust height
plot_resolution <- 300    # DPI for high resolution

# Save plot as PNG for Feature Importance
ggsave(
  filename = here("Q:/Dante/ZoopMetaB/rfr_attempt/rfr_coi_otu_importance_plot.png"),
  plot = ggplot(best_otu_importance, aes(x = reorder(Feature_label, Importance), y = Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste("Feature Importance for \n Best Model (", best_result$n_otus, "OTUs)", sep = ""),
      x = "OTUs/Species",
      y = "Importance"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 10, family = "serif"),  # Set font size and style
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold")  # Center title
    ),
  width = plot_width,
  height = plot_height,
  dpi = plot_resolution
)



# 18S ---------------------------------------------------------------------

####### 18S

#18s reads
zhan_otu=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash")%>%
  select(where(~ !is.na(.[[1]])))
rownames(zhan_otu) <- paste0("otu_", gsub(" ", ".", rownames(zhan_otu)))

zhan_meta=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  column_to_rownames("Sample_ID_dot") %>%
  dplyr::select(-X)
zhan_taxa=read.csv(here("data/phyloseq_bio_data/18s/metazoopruned18s_tax.csv")) %>% column_to_rownames("Hash") %>%
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family))
rownames(zhan_taxa) <- paste0("otu_", gsub(" ", ".", rownames(zhan_taxa)))



#Convert to phyloseq

OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(zhan_meta)
Phy_merged_18s <- phyloseq(OTU, TAX, meta)

# Extract OTU and taxa data
otu_data <- as(otu_table(Phy_merged_18s), "matrix")
taxa_data <- as.data.frame(tax_table(Phy_merged_18s))

# Replace spaces in taxa table columns
taxa_data <- taxa_data %>%
  mutate(across(everything(), ~ gsub(" ", ".", .)))  # Replace spaces with "."
meta_data <- as(sample_data(Phy_merged_18s), "data.frame") # Metadata

# # Step 3: Use species names as OTU table rownames
# species_names <- taxa_data$
# rownames(otu_data) <- species_names

# Normalize OTU table (relative abundance)
otu_normalized <- apply(otu_data, 2, function(x) x / sum(x))

# Step 4: Align with metadata and prepare for Random Forest
otu_samples <- colnames(otu_normalized)
meta_samples <- rownames(meta_data)
common_samples <- intersect(otu_samples, meta_samples)

otu_normalized <- otu_normalized[, common_samples]
meta_data <- meta_data[common_samples, ]

# Check for and handle missing values in PC1
if (any(is.na(meta_data$PC1))) {
  meta_data <- meta_data[!is.na(meta_data$PC1), ]
}

# Prepare data frame for modeling
rownames(otu_normalized) <- gsub(" ", ".", rownames(otu_normalized))
otu_df <- as.data.frame(t(otu_normalized)) # Transpose for predictors

# Ensure unique column names
colnames(otu_df) <- make.unique(colnames(otu_df))

# Merge columns containing "Species" using dplyr
otu_df <- otu_df %>%
  mutate(unidentified.Family = rowSums(select(., contains("Family"))))%>%
  select(-contains("Family"))

otu_df$PC1 <- meta_data$PC1

# Step 5: Train-Test Split
set.seed(123)
train_index <- createDataPartition(otu_df$PC1, p = 0.8, list = FALSE)
train_data <- otu_df[train_index, ]
test_data <- otu_df[-train_index, ]

# Step 6: Fit Random Forest Regression
set.seed(123)
rf_model <- randomForest(
  PC1 ~ ., 
  data = train_data,
  ntree = 500, 
  importance = TRUE
)

# Step 7: Evaluate Model
# Extract feature columns from training data
train_features <- train_data %>% select(-PC1)
test_features <- test_data %>% select(-PC1)

# Ensure columns match between training and test sets
if (!all(colnames(train_features) == colnames(test_features))) {
  stop("Column names in test data do not match the training data!")
}

# Predict on test data
predicted <- predict(rf_model, newdata = test_features)

# Compare predictions with actual values
actual <- test_data$PC1

# Calculate RMSE
rmse <- sqrt(mean((predicted - actual)^2))
cat("Root Mean Squared Error (RMSE):", rmse, "\n")

# Plot model performance
performance_plot <- data.frame(Actual = actual, Predicted = predicted)
ggplot(performance_plot, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Model Performance: Predicted vs Actual",
    x = "Actual PC1",
    y = "Predicted PC1"
  ) +
  theme_minimal()

# Step 8: Feature Importance
importance <- importance(rf_model)
importance_df <- data.frame(Feature = rownames(importance), Importance = importance[, 1]) %>%
  arrange(desc(Importance))

# Join importance_df with taxa_data to get species names
importance_with_taxa <- importance_df %>%
  left_join(taxa_data%>% rownames_to_column(var = "Feature"), by = c("Feature")) # Ensure taxa_data has row names or a matching column

# Replace spaces with "." in species names before handling duplicates
importance_with_taxa <- importance_with_taxa %>%
  mutate(Family = gsub(" ", ".", Family),
         Family = gsub(" ", ".", Family))%>%
  # Replace Feature column with Family names, or use Family names if Family is missing
  mutate(Feature_label = ifelse(!is.na(Family), paste0(Family,"_",Feature), paste0("unidentified_",Order,"_",Feature)))





# Filter out "unidentified_NA_NA" features
importance_with_taxa <- importance_with_taxa %>%
  filter(Feature != "unidentified_NA_NA")

# Order Feature by Importance for plotting
importance_with_taxa <- importance_with_taxa %>%
  mutate(Feature_label = factor(Feature_label, levels = Feature_label[order(Importance)]))

# Plot top 20 important features with species names
top_features <- head(importance_with_taxa, 20)
ggplot(top_features, aes(x = Feature_label, y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 Important Species Predicting PC1", x = "Species", y = "Importance") +
  theme_minimal()

# Optional: Extract the core species names or OTUs
core_otus <- top_features$Feature
cat("Core OTUs/Species:", paste(core_otus, collapse = ", "), "\n")



# Step 9: Loop through subsets of OTUs to find the best model
# Initialize storage for results
subset_results <- list()

# Rank OTUs by importance
ranked_otus <- importance_with_taxa$Feature

# Loop through decreasing numbers of OTUs
for (n_otus in seq(length(ranked_otus), 10, by = -10)) {
  
  # Subset top n OTUs
  selected_otus <- head(ranked_otus, n_otus)
  
  # Filter training and test datasets to include only the selected OTUs
  train_subset <- train_data %>% select(all_of(selected_otus), PC1)
  test_subset <- test_data %>% select(all_of(selected_otus), PC1)
  
  # Fit Random Forest model on the subset
  set.seed(123)
  rf_subset_model <- randomForest(
    PC1 ~ ., 
    data = train_subset,
    ntree = 500
  )
  
  # Predict and calculate RMSE
  subset_predicted <- predict(rf_subset_model, newdata = test_subset %>% select(-PC1))
  subset_actual <- test_subset$PC1
  subset_rmse <- sqrt(mean((subset_predicted - subset_actual)^2))
  
  
  # Store results
  subset_results[[as.character(n_otus)]] <- list(
    n_otus = n_otus,
    model = rf_subset_model,
    rmse = subset_rmse
  )
}

# Step 10: Summarize results and find the optimal model
results_df <- do.call(rbind, lapply(subset_results, function(res) {
  data.frame(n_otus = res$n_otus, RMSE = res$rmse)
}))

# Find the best model (min RMSE)
best_result <- results_df[which.min(results_df$RMSE), ]
cat("Best model uses", best_result$n_otus, "OTUs with RMSE:", best_result$RMSE, "\n")

# Plot RMSE vs. number of OTUs
ggplot(results_df, aes(x = n_otus, y = RMSE)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Model Performance vs Number of OTUs",
    x = "Number of OTUs",
    y = "RMSE"
  ) +
  theme_minimal()



# Find the best result (minimum RMSE)
best_result <- subset_results[[as.character(best_result$n_otus)]]

# Extract the best model and data
best_model <- best_result$model
best_test_subset <- test_data %>% select(all_of(head(ranked_otus, best_result$n_otus)), PC1)

# Predict using the best model
best_predicted <- predict(best_model, newdata = best_test_subset %>% select(-PC1))
best_actual <- best_test_subset$PC1

# Remake model performance plot
best_performance_plot <- data.frame(Actual = best_actual, Predicted = best_predicted)
ggplot(best_performance_plot, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(
    title = paste("Model Performance with Best Result (", best_result$n_otus, "OTUs)", sep = ""),
    x = "Actual PC1",
    y = "Predicted PC1"
  ) +
  theme_minimal()



# Evaluate the best model objectively
subset_rmse <- sqrt(mean((best_predicted - best_actual)^2))  # RMSE
subset_mae <- mean(abs(best_predicted - best_actual))        # MAE
subset_r2 <- 1 - (sum((best_actual - best_predicted)^2) / sum((best_actual - mean(best_actual))^2))  # R-squared

cat("Model Evaluation Metrics:\n")
cat("RMSE:", subset_rmse, "\n")
cat("MAE:", subset_mae, "\n")
cat("R-squared:", subset_r2, "\n")


# Get the top OTUs contributing to the best model
best_contributing_otus <- head(ranked_otus, best_result$n_otus)

# Filter importance_with_taxa for the top contributing OTUs
best_otu_importance <- importance_with_taxa %>%
  filter(Feature %in% best_contributing_otus)%>%
  mutate(Feature_label = substr(Feature_label, 1, 40))  %>%  # Truncate Feature to 40 characters
  arrange(desc(Importance)) %>%  # Sort by Importance in descending order
  slice_head(n = 30)  # Select the top 30 OTUs

# Plot the importance of contributing OTUs
ggplot(best_otu_importance, aes(x = reorder(Feature_label, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = paste("Feature Importance for Best Model (", best_result$n_otus, "OTUs)", sep = ""),
    x = "OTUs/Species",
    y = "Importance"
  ) +
  theme_minimal()


# Define plot dimensions and resolution for L&O
plot_width <- 12  # Convert cm to inches (single-column width)
plot_height <- 4  # Adjust height
plot_resolution <- 300    # DPI for high resolution

# Save plot as PNG for Feature Importance
ggsave(
  filename = here("Q:/Dante/ZoopMetaB/rfr_attempt/rfr_18s_otu_importance_plot.png"),
  plot = ggplot(best_otu_importance, aes(x = reorder(Feature_label, Importance), y = Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste("Feature Importance for \n Best Model (", best_result$n_otus, "OTUs)", sep = ""),
      x = "OTUs/Species",
      y = "Importance"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 10, family = "serif"),  # Set font size and style
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold")  # Center title
    ),
  width = plot_width,
  height = plot_height,
  dpi = plot_resolution
)


