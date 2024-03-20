#Exploratory environmental data analysis
library(phyloseq)
library(tidyverse)
library(fido)
library(here)

here()
#COI reads
otucoi=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash")%>%
  dplyr::select(where(~ !is.na(.[[1]]))) %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)
taxcoi=read.csv(here("data/phyloseq_bio_data/COI/coi_taxa_table_eDNA_metazoogene.csv")) %>%
  column_to_rownames("X")
taxcoi=tax_table(as.matrix(taxcoi))
metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size))

dat_all=phyloseq(otucoi,taxcoi,metacoi)
dat=merge_samples(dat_all,"Sample_ID_short",fun= mean)%>%
  filter_taxa(function(x) sum(x > 3) > 0.10*length(x), TRUE)

set.seed(899)


##Env Correlation matrix
meta_corr=metacoi %>% dplyr::select(-Sample_ID_short) %>%
  cor(.)

# Compute p-values using correlation matrix
p_values <- cor.mtest(metacoi %>% dplyr::select(-Sample_ID_short))$p %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(cols = -variable, names_to = "variable2", values_to = "p.value")%>%
  # Adjust p-values using Benjamini-Hochberg correction
  mutate(p.adj= p.adjust(p.value, method = "BH") )

p_vals_adj=p_values %>%
  group_by(variable, variable2) %>%
  summarise(p.adj = mean(p.adj, na.rm = TRUE)) %>%
  ungroup() %>% # Ensure to ungroup the data after summarizing 
  pivot_wider(names_from = variable, values_from = p.adj) %>%
  column_to_rownames("variable2") %>%
  as.matrix()




# Corr using corrplot --------------------------------------------------------------------

# Sort the row names and column names to ensure alignment
meta_corr <- meta_corr[order(rownames(meta_corr)), order(colnames(meta_corr))]
p_vals_adj <- p_vals_adj[order(rownames(p_vals_adj)), order(colnames(p_vals_adj))]


corr_plot=corrplot(meta_corr,p.mat=p_vals_adj, type = 'lower', order = 'FPC', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45)

#PNG & PDF Save
ggsave(
  filename = here("plots/Q1_physical_analysis/corr_plot_p_adj.png"),
  plot = corr_plot,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/Q1_physical_analysis/corr_plot_p_adj.pdf"),
  plot = corr_plot,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


