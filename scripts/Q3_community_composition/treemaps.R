#Script for Question 3: Community Composition Using Treemaps for RRA


#Load packages
library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
library(ggpubr)
library(fido)
library(phyloseq)
here()

#Load functions
source(("scripts/helpful_functions/treemap_funs_Capone.R"))


############## COI


#COI reads
leray_metazoo_otucoi=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash")%>%
  select(where(~ !is.na(.[[1]])))

leray_metazoo_meta=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  column_to_rownames("Sample_ID_dot") %>%
  dplyr::select(-X)
leray_metazoo_taxa=read.csv(here("data/phyloseq_bio_data/COI/coi_taxa_table_eDNA_metazoogene.csv")) %>% column_to_rownames("X")


#Convert to phyloseq

OTU = otu_table(as.matrix(leray_metazoo_otucoi), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(leray_metazoo_taxa))
meta=sample_data(leray_metazoo_meta)
Phy_merged_coi <- phyloseq(OTU, TAX, meta)

#### Transform to Long
phy_merged_long_coi=phyloseq_transform_to_long((Phy_merged_coi)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species") 

#### Transform to Long
phy_norm_coi=phyloseq_transform_to_long(phyloseq_normalize_median(Phy_merged_coi)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", NA, Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family))


## ALL TOP 15
p_coi=phyloseq_long_treemap_top15(phy_norm_coi, Family, Genus ,"COI All",colors=NULL, label_group1 = TRUE)
p_coi

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/coi_top15.png"),
  plot = p_coi,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/coi_top15.pdf"),
  plot = p_coi,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

###By onshore offshore
off_on=unique(phy_merged_long_coi$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_merged_long_coi[phy_merged_long$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Species,Genus,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_coi_on=list.plots[[1]]
p_coi_on

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_coi_top10.png"),
  plot = p_coi_on,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_coi_top10.pdf"),
  plot = p_coi_on,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_coi_off=list.plots[[2]]
p_coi_off

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_coi_top10.png"),
  plot = p_coi_off,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_coi_top10.pdf"),
  plot = p_coi_off,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)







####### 18S

#COI reads
zhan_otu=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash")%>%
  select(where(~ !is.na(.[[1]])))

zhan_meta=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  column_to_rownames("Sample_ID_dot") %>%
  dplyr::select(-X)
zhan_taxa=read.csv(here("data/phyloseq_bio_data/18s/metazoopruned18s_tax.csv")) %>% column_to_rownames("Hash") %>%
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family))


#Convert to phyloseq

OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(zhan_meta)
Phy_merged_18s <- phyloseq(OTU, TAX, meta)

#### Transform to Long
phy_merged_long_18s=phyloseq_transform_to_long((Phy_merged_18s)) %>%
  filter(Order != "Order")%>%
  filter(Family != "Family") %>%
  filter(Genus != "Genus") %>%
  mutate(Species = ifelse(Species == "", NA, Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family))

#### Transform to Long

phy_norm_18s=phyloseq_transform_to_long((Phy_merged_18s)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")
group1=phy_norm$Genus
group2=phy_norm$Species


# All
p_18s=phyloseq_long_treemap_top(phy_merged_long_18s, Family,Genus,"18S All",top=15, label_group1 = TRUE)
p_18s
#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/zhan_top15.png"),
  plot = p_18s,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/zhan_top15.pdf"),
  plot = p_18s,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

###By onshore offshore
off_on=unique(phy_merged_long_18s$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_merged_long_18s[phy_merged_long$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Genus,Family,titles[i],top=10,colors=NULL, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_18s_on=list.plots[[1]]
p_18s_on

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_18s_top10.png"),
  plot = p_18s_on,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_18s_top10.pdf"),
  plot = p_18s_on,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_18s_off=list.plots[[2]]
p_18s_off

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_18s_top10.png"),
  plot = p_18s_off,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_18s_top10.pdf"),
  plot = p_18s_off,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)



###============ Stacked Bar plot
phy_merged_long_18s_plot=phy_merged_long_18s%>%
  group_by(cycle) %>%
  mutate(prop = n_reads / sum(n_reads))%>%
  group_by(Family,cycle)
  filter(prop >= 1)

# Identify families with proportion < 5%
low_prop_families <- phy_merged_long_18s %>%
  group_by(cycle) %>%
  filter(prop < 0.001) %>%
  group_by(cycle)  %>%
  summarise(across(c(n_reads, prop), sum),
            across(where(is.numeric) & !all_of(c("n_reads", "prop")), mean),
            across(where(is.character), first)) %>%
  mutate(Family="Other")




phy_merged_long_18s_plot %>%
  # bind_rows(.,low_prop_families) %>%
  # filter(prop > 0.01 | Family=="Other") %>% # Create a new data frame with 'other' category
  ggplot(aes(x = as.factor(cycle), y = prop, fill = Family)) +
  geom_bar(stat = "identity") +
  labs(x = "Cycle", y = "Proportion of Total Reads", fill = "Family") +
  ggtitle("Stacked Barplot of Proportions by Family and Cycle") +
  theme_minimal() 
