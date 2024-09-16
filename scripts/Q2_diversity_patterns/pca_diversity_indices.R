#PCA Script

#Load packages and set path
library(tidyverse)
library(here)
library(lubridate)
library(dplyr)
library(matrixStats)
library(ggpubr)
library(stringr)
library(phyloseq)
library(tibble)
library(tidyr)
library(gridExtra)
library(extrafont)
here("")
#Add functions for myself
source(("Zoop_Patterns/scripts/helpful_functions/treemap_funs_Capone.R"))
source("Zoop_Patterns/scripts/helpful_functions/phyloseq_mapping_funs.R")
source("Zoop_Patterns/scripts/helpful_functions/general_helper_functions.R")



#Part one: Read in the physical environmental data

env_metadata<-read.csv(here("Zoop_Patterns/data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv")) %>% dplyr::select(-c("X"))%>%
  column_to_rownames("Sample_ID_dot") %>%
  mutate(PC1=PC1*-1)

#Load in the phyloseq data and format to a table

#COI raw reads
coi_metazoo_otu=read.csv(here("Zoop_Patterns/data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
coi_metazoo_meta=env_metadata 
  
coi_metazoo_taxa=read.csv(here("Zoop_Patterns/data/phyloseq_bio_data/COI/metazooprunedcoi_tax.csv")) %>% column_to_rownames("Hash")


#Merge by site
coi_metazoo_meta_all=coi_metazoo_meta %>% group_by(Sample_ID_short) %>%
  summarize_all(median) %>% 
  column_to_rownames("Sample_ID_short")

OTU = otu_table(as.matrix(coi_metazoo_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_metazoo_taxa))
meta=sample_data(coi_metazoo_meta)
meta$cycle=meta$cycle %>% as.factor()
Phy_coi_raw <- phyloseq(OTU, TAX, meta)


#Merge by sample site
Phy_merged_coi_raw <- merge_samples(Phy_coi_raw,c("Sample_ID_short"))


#Issues with metadata...create a new phyloseq object with merged metadata
Phy_merged_coi_raw=phyloseq(otu_table(Phy_merged_coi_raw),tax_table(Phy_merged_coi_raw),sample_data(coi_metazoo_meta_all))

#18s
zhan_metazoo_otu=read.csv(here("Zoop_Patterns/data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
zhan_metazoo_meta=env_metadata
zhan_metazoo_taxa=read.csv(here("Zoop_Patterns/data/phyloseq_bio_data/18S/metazoopruned18s_tax.csv")) %>% column_to_rownames("Hash")


#Merge by site
zhan_metazoo_meta_all=zhan_metazoo_meta %>% group_by(Sample_ID_short) %>%
  summarize_all(median) %>% 
  column_to_rownames("Sample_ID_short")

OTU = otu_table(as.matrix(zhan_metazoo_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_metazoo_taxa))
meta=sample_data(zhan_metazoo_meta)
meta$cycle=meta$cycle %>% as.factor()
Phy_zhan_raw <- phyloseq(OTU, TAX, meta)

#Merge by sample site
Phy_merged_zhan_raw <- merge_samples(Phy_zhan_raw,c("Sample_ID_short"))


#Issues with metadata...create a new phyloseq object with merged metadata
Phy_merged_zhan_raw=phyloseq(otu_table(Phy_merged_zhan_raw),tax_table(Phy_merged_zhan_raw),sample_data(zhan_metazoo_meta_all))





###Correlate PC1 and Shannon

#compute Shannon & Chao index
plot_data_coi=coi_metazoo_meta_all %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID_short") %>%
  mutate(estimate_richness(Phy_merged_coi_raw, measures="Shannon")) %>% 
  mutate(estimate_richness(Phy_merged_coi_raw, measures="Chao1"))

plot_data_18s=zhan_metazoo_meta_all %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID_short") %>%
  mutate(estimate_richness(Phy_merged_zhan_raw, measures="Shannon")) %>% 
  mutate(estimate_richness(Phy_merged_zhan_raw, measures="Chao1"))



#Correlate
#Test corr
# Calculate Pearson's correlation coefficient and p-value
correlation_result_coi <- cor.test(plot_data_coi$PC1, plot_data_coi$Shannon, method = "spearman")
correlation_result_18s <- cor.test(plot_data_18s$PC1, plot_data_18s$Shannon, method = "pearson")

correlation_result_coi
correlation_result_18s
# Extract Pearson's r and p-value
pearsons_r <- correlation_result_coi$estimate
p_value <- correlation_result_coi$p.value

pearsons_r <- correlation_result_18s$estimate
p_value <- correlation_result_18s$p.value

# Print Pearson's r and p-value
print(paste("Pearson's r:", round(pearsons_r, 3)))
print(paste("p-value:", format(p_value, scientific = FALSE)))


# Run linear regression
lm_model <- lm(Shannon ~ PC1, data = plot_data_18s)

# Summary of linear regression
lm_summary <- summary(lm_model)

# Extracting R-squared and p-value
r_squared <- lm_summary
p_value <- lm_summary$coefficients[2, 4]



#Correlation calculations
lm_model_coi <- lm(Shannon ~ PC1, data = plot_data_coi)

# Summary of linear regression
summary(lm_model_coi)

## Create a scatter plot with regression line, confidence intervals, and color by 'cycle'
#My colors for Cycles
txt_sz=24
my_palette=custom_pallete()
coi_plot=ggplot(plot_data_coi, aes(x = PC1, y = Shannon)) +
  geom_point(size=8, aes(shape=cycle, fill=cycle))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  coord_cartesian(ylim = c(1.5, 4), xlim = c(min(plot_data_18s$PC1), 7))+
  scale_x_continuous(breaks = seq(-6, 7, by = 2))+
  scale_fill_manual(values = my_palette) +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="COI",
       shape="Cycle", fill="Cycle") +
  scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  stat_cor(method="pearson", label.x = 0, label.y = 3.5, size=8)+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 4, label.y = 3.7)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))
coi_plot
saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_all.pdf"), 
    plot = coi_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_all.png"), 
    plot = coi_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

#18s all
zhan_plot=ggplot(plot_data_18s, aes(x = PC1, y = Shannon)) +
  geom_point(size=8, aes(shape=cycle, fill=cycle))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  coord_cartesian(ylim = c(0, 4), xlim = c(min(plot_data_18s$PC1), 7))+
  scale_x_continuous(breaks = seq(-6, 7, by = 2))+
  scale_fill_manual(values = my_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))+
  # scale_y_continuous(breaks = seq(0, 4, length.out = 10))+
  scale_fill_manual(values = my_palette) +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="18S",
       shape="Cycle", fill="Cycle") +
  scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  stat_cor(method="pearson", label.x = 0, label.y = 3.5, size=8)+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 4, label.y = 3.7)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))
zhan_plot
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_all.pdf"), 
    plot = zhan_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_all.png"), 
    plot = zhan_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }





#Both
both_plot=grid.arrange(coi_plot,zhan_plot)
both_plot

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all.pdf"), 
    plot = both_plot,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all.png"), 
    plot = both_plot,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }


#Save
#COI
output_path <- here("plots", "Q1_physical_analysis")
ggsave(file.path(output_path, "pc1_vs_shannon_coi.png"), coi_plot, width = 10, height = 6, units = "in")

#18s
output_path <- here("plots", "Q1_physical_analysis")
ggsave(file.path(output_path, "pc1_vs_shannon_18s.png"), zhan_plot, width = 10, height = 6, units = "in")

#Both plots 
output_path <- here("plots", "Q1_physical_analysis")
ggsave(file.path(output_path, "pc1_vs_shannon_both.png"), both_plot, width = 10, height = 6, units = "in")






#### PC1 vs Shannon for each size
#compute Shannon index
shannon_coi=estimate_richness(Phy_coi_raw, measures="Shannon") %>% 
  rownames_to_column("Sample_ID")
shannon_18s=estimate_richness(Phy_zhan_raw, measures="Shannon") %>% 
  rownames_to_column("Sample_ID")


plot_data_coi=coi_metazoo_meta %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(.,shannon_coi, by="Sample_ID")%>%
  mutate(max_size=as.factor(max_size))

plot_data_18s=zhan_metazoo_meta %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(.,shannon_18s, by="Sample_ID")%>%
  mutate(max_size=as.factor(max_size))

#Test corr
# Calculate Pearson's correlation coefficient and p-value
correlation_result <- cor.test(plot_data_coi$PC1, plot_data_coi$Shannon, method = "pearson")

# Extract Pearson's r and p-value
pearsons_r <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print Pearson's r and p-value
print(paste("Pearson's r:", round(pearsons_r, 3)))
print(paste("p-value:", format(p_value, scientific = TRUE)))

#
facet_correlation <- plot_data_ %>%
  group_by(max_size) %>%
  summarise(pearsons_r = cor(PC1, Shannon, method = "pearson"), 
            p_value = cor.test(PC1, Shannon, method = "pearson")$p.value)

# Print the result
print(facet_correlation)

txt_sz=24
coi_plot_sized=ggplot(plot_data_coi, aes(x = PC1, y = Shannon)) +
  geom_point(size = 6, aes(shape = cycle, fill = cycle), show.legend = TRUE) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, aes(color = max_size), show.legend = FALSE, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 5), xlim = c(min(plot_data_coi$PC1),  max(plot_data_coi$PC1))) +
  scale_x_continuous(breaks = seq(-6, max(plot_data_coi$PC1), by = 1.5)) +
  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="COI") +
  # scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  stat_cor(method="pearson", label.x = 2, label.y = 4,
           size=8)+
  theme_classic()+
  facet_wrap(~max_size, nrow = 3, labeller = labeller(max_size = c("0.5" = "0.2-0.5 mm", "1" = "0.5-1 mm", "2" = "1-2 mm"))) +
  # Altering font sizes
  theme(strip.text = element_text(size = txt_sz),
        axis.text.x = element_text(size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title.x = element_text(size = txt_sz),
        axis.title.y = element_text(size = txt_sz)) +
  guides(color = "none") 
coi_plot_sized

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_sized.pdf"), 
    plot = coi_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_sized.png"), 
    plot = coi_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }


#Test corr
# Calculate Pearson's correlation coefficient and p-value
correlation_result <- cor.test(plot_data_18s$PC1, plot_data_18s$Shannon, method = "pearson")

# Extract Pearson's r and p-value
pearsons_r <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print Pearson's r and p-value
print(paste("Pearson's r:", round(pearsons_r, 3)))
print(paste("p-value:", format(p_value, scientific = TRUE)))

#
facet_correlation <- plot_data_18s %>%
  group_by(max_size) %>%
  summarise(pearsons_r = cor(PC1, Shannon, method = "pearson"), 
            p_value = cor.test(PC1, Shannon, method = "pearson")$p.value)

# Print the result
print(facet_correlation)


#Plot to visualize
zhan_plot_sized <- ggplot(plot_data_18s, aes(x = PC1, y = Shannon)) +
  geom_point(size = 6, aes(shape = cycle, fill = cycle), show.legend = TRUE) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, aes(color = max_size), show.legend = FALSE, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 4), xlim = c(min(plot_data_18s$PC1), max(plot_data_18s$PC1)+1)) +
  scale_x_continuous(breaks = seq(-6, max(plot_data_18s$PC1)+1, by = 1.5)) +
  scale_y_continuous(breaks = seq(0, 4, by = 1)) +
  scale_fill_manual(values = my_palette) +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="18S") +
  # scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  stat_cor(method="pearson", label.x = 4, label.y = 3.5,
  size=8)+
  theme_classic()+
  facet_wrap(~max_size, nrow = 3, labeller = labeller(max_size = c("0.5" = "0.2-0.5 mm", "1" = "0.5-1 mm", "2" = "1-2 mm"))) + 
  # Altering font sizes
  theme(strip.text = element_text(size = txt_sz),
        axis.text.x = element_text(size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title.x = element_text(size = txt_sz),
        axis.title.y = element_text(size = txt_sz)) +
  guides(color = "none") 
zhan_plot_sized

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_sized.pdf"), 
    plot = zhan_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_sized.png"), 
    plot = zhan_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }


both_plot_sized=grid.arrange(coi_plot_sized,zhan_plot_sized, nrow=2)
both_plot_sized

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all_sized.pdf"), 
    plot = both_plot_sized,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all_sized.png"), 
    plot = both_plot,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }

