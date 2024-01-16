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
here()

#Part one: Read in the physical environmental data

env_metadata<-read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv")) %>% dplyr::select(-c("X"))%>%
  column_to_rownames("Sample_ID_dot")

#Load in the phyloseq data and format to a table

#COI raw reads
coi_metazoo_otu=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
coi_metazoo_meta=env_metadata
coi_metazoo_taxa=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_tax.csv")) %>% column_to_rownames("Hash")


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
zhan_metazoo_otu=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
zhan_metazoo_meta=env_metadata
zhan_metazoo_taxa=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_tax.csv")) %>% column_to_rownames("Hash")


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

#compute Shannon index
shannon_coi=estimate_richness(Phy_merged_coi_raw, measures="Shannon") %>% 
  rownames_to_column("Sample_ID_short")
shannon_18s=estimate_richness(Phy_merged_zhan_raw, measures="Shannon") %>% 
  rownames_to_column("Sample_ID_short")

#Chao
chao_coi=estimate_richness(Phy_merged_coi_raw, measures="Chao1") %>% 
  rownames_to_column("Sample_ID_short")
chao_zhan=estimate_richness(Phy_merged_zhan_raw, measures="Chao1") %>% 
  rownames_to_column("Sample_ID_short")


plot_data_coi=coi_metazoo_meta_all %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID_short") %>%
  left_join(.,shannon_coi, by="Sample_ID_short")

plot_data_18s=zhan_metazoo_meta_all %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID_short") %>%
  left_join(.,shannon_18s, by="Sample_ID_short")

#Correlate


## Create a scatter plot with regression line, confidence intervals, and color by 'cycle'

coi_plot=ggplot(plot_data_coi, aes(x = PC1, y = Shannon, color=cycle)) +
  geom_point(size=6, aes(shape=cycle)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  labs(x = "PC1", y = "Shannon Diversity Index") +
  scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  theme_minimal()

zhan_plot=ggplot(plot_data_18s, aes(x = PC1, y = Shannon, color=cycle)) +
  geom_point(size=6, aes(shape=cycle)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  labs(x = "PC1", y = "Shannon Diversity Index") +
  scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  theme_minimal()


grid.arrange(coi_plot,zhan_plot)
