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


# Convert the correlation matrix into a long format
cor_long <- as.data.frame(as.table(meta_corr))

# Plot using ggplot2
ggplot(data = cor_long, aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=Freq), color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white", 
                       midpoint=0, limit=c(-1,1), space="Lab", 
                       name="Correlation") +
  geom_text(aes(label=sprintf("%.2f", Freq)), vjust=1) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust=1))


## Violin plots


