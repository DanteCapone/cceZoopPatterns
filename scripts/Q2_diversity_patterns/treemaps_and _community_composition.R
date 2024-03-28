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

# COI ---------------------------------------------------------------------



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
phy_coi_majority=phyloseq_transform_to_long((Phy_merged_coi)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", NA, Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family)) %>%
  filter((Order %in% c("Calanoida","Euphausiacea")))

phy_coi_minority=phyloseq_transform_to_long((Phy_merged_coi)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", NA, Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family)) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea")))


## ALL Majority taxa
p_coi_maj=phyloseq_long_treemap_top(phy_coi_majority, Genus,Family ,"COI All",20,colors=NULL, label_group1 = TRUE)
p_coi_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/coi_majority.png"),
  plot = p_coi_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/coi_majority.pdf"),
  plot = p_coi_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


## ALL Minority taxa
p_coi_min=phyloseq_long_treemap_top(phy_coi_minority, Family, Genus ,"COI All",20,colors=NULL, label_group1 = TRUE)
p_coi_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/coi_minority.png"),
  plot = p_coi_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/coi_majority.pdf"),
  plot = p_coi_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)





###By onshore offshore minority coi
off_on=unique(phy_merged_long_coi$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_coi_minority[phy_coi_minority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Species,Genus,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_coi_on_min=list.plots[[1]]
p_coi_on_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_coi_minority.png"),
  plot = p_coi_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_coi_minority.pdf"),
  plot = p_coi_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_coi_off_min=list.plots[[2]]
p_coi_off_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_coi_minority.png"),
  plot = p_coi_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_coi_minority.pdf"),
  plot = p_coi_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


###By onshore offshore majority coi
off_on=unique(phy_merged_long_coi$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_coi_majority[phy_coi_majority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Species,Genus,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_coi_on_maj=list.plots[[1]]
p_coi_on_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_coi_majority.png"),
  plot = p_coi_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_coi_majority.pdf"),
  plot = p_coi_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_coi_off_maj=list.plots[[2]]
p_coi_off_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_coi_majority.png"),
  plot = p_coi_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_coi_majority.pdf"),
  plot = p_coi_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)






# 18s ---------------------------------------------------------------------


####### 18S

#18s reads
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
phy_18s_majority=phyloseq_transform_to_long((Phy_merged_18s)) %>%
  mutate(Genus = ifelse(is.na(Genus), "unknown Genus", Genus)) %>%
  mutate(Genus = ifelse(Genus == "", "unknown Genus", Genus))%>%
  mutate(Family = ifelse(is.na(Family), "unknown Family", Family)) %>%
  filter((Order %in% c("Calanoida","Euphausiacea"))) 
  

phy_18s_minority=phyloseq_transform_to_long((Phy_merged_18s)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", "NA", Species))%>%
  mutate(Genus = ifelse(is.na(Genus), "unknown Genus", Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family)) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea")))


## ALL Majority taxa
p_18s_maj=phyloseq_long_treemap_top(phy_18s_majority, Genus,Family ,"18s All",20,colors=NULL, label_group1 = TRUE)
p_18s_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/18s_majority.png"),
  plot = p_18s_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/18s_majority.pdf"),
  plot = p_18s_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


## ALL Minority taxa
p_18s_min=phyloseq_long_treemap_top(phy_18s_minority, Genus,Family ,"18s All",20,colors=NULL, label_group1 = TRUE)
p_18s_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/18s_minority.png"),
  plot = p_18s_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/18s_majority.pdf"),
  plot = p_18s_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)







###By onshore offshore minortiy
off_on=unique(phy_18s_majority$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_18s_minority[phy_18s_minority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Genus,Family,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_18s_on_min=list.plots[[1]]
p_18s_on_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_18s_minority_fam.png"),
  plot = p_18s_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_18s_minority_fam.pdf"),
  plot = p_18s_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_18s_off_min=list.plots[[2]]
p_18s_off_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_18s_minority_fam.png"),
  plot = p_18s_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_18s_minority_fam.pdf"),
  plot = p_18s_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


###By onshore offshore majority
off_on=unique(phy_18s_majority$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_18s_majority[phy_18s_majority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Genus,Family,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_18s_on_maj=list.plots[[1]]
p_18s_on_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_18s_majority_fam.png"),
  plot = p_18s_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_18s_majority_fam.pdf"),
  plot = p_18s_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_18s_off_maj=list.plots[[2]]
p_18s_off_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_18s_majority_fam.png"),
  plot = p_18s_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_18s_majority_fam.pdf"),
  plot = p_18s_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)



# Stacked Bar Plots -------------------------------------------------------



# COI ---------------------------------------------------------------------

Phy_glom_coi <- Phy_merged_coi %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop))

Phy_glom_coi_minority <- Phy_merged_coi %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop)) %>% 
  ungroup()





#Reorder Calalnoida and Euphausiacea for plotting
Phy_glom_coi %>%
  ungroup() %>% 
  select(Order) %>% 
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>% 
  unique() %>% 
  as.matrix()->other_orders 

Phy_glom_coi=Phy_glom_coi%>%
  mutate(Order=factor(Order, levels =c("Calanoida","Euphausiacea",other_orders)))

#PC1 Labels
labels_for_PC1=Phy_glom_coi %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

#Taxa color maps
orders_coi=unique(Phy_glom_coi$Order)

# Define two Brewer palettes
palette1 <- brewer.pal(8, "Accent")
palette2 <- brewer.pal(8, "Pastel1")

# Join the two palettes
contrast_palette <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#9edae5" , "#9467bd", "#8c564b", 
                       "#17becf", "#7f7f7f", "#bcbd22", "#e377c2", "#aec7e8", "#ffbb78",
                       "#98df8a", "#ff9896", "#ffaec9","#fc798a")
palette_named <- setNames(contrast_palette, orders_coi)


Phy_glom_coi %>%
  group_by(Order) %>%
  # filter(Order %in% c("Calanoida","Euphausiacea")) %>%
  # filter(prop < 0.001)  %>%
  # bind_rows(.,low_prop_families) %>%
  # filter(prop > 0.01 | Family=="Other") %>% # Create a new data frame with 'other' category
  ggplot(.,aes(x = as.factor(PC1), y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Proportion of Total Reads", fill = "Order") +
  # ggtitle("Raw Relative Read Abundances By Cycle and Family") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_named) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)-> majority_p 



majority_p

Phy_glom_coi %>%
  group_by(Order) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  # filter(prop < 0.001)  %>%
  # bind_rows(.,low_prop_families) %>%
  # filter(prop > 0.01 | Family=="Other") %>% # Create a new data frame with 'other' category
  ggplot(.,aes(x = as.factor(PC1), y = asin(sqrt((prop))), fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Arcsine-Squareroot\nProportion of Residual Reads", fill = "Order") +
  # ggtitle("Raw Relative Read Abundances By Cycle and Family") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_named) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) -> minority_p 
minority_p
stacked_bar_coi=gridExtra::grid.arrange(majority_p,minority_p,ncol=1)


#PNG & PDF Save
ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_coi.png"),
  plot = stacked_bar_coi,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)

ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_coi.pdf"),
  plot = stacked_bar_coi,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)




# 18S ---------------------------------------------------------------------


Phy_glom_18s <- Phy_merged_18s %>%
  tax_glom(taxrank="Family") %>%
  phyloseq_transform_to_long(.) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop))

Phy_glom_18s_minority <- Phy_merged_18s %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop))



# Match categories and get colors for overlapping categories
matching_categories <- intersect(orders_coi,orders_18s)
matching_colors <- palette_named [match(matching_categories, orders_coi)]

# New categories in the second dataset
new_categories <- setdiff(orders_18s, matching_categories)


# Generate new colors for new categories
new_colors <- brewer.pal(length(new_categories), "Pastel2")

# Combine colors for both datasets
palette_18s <- c(setNames(matching_colors, matching_categories), setNames(new_colors, new_categories))



#Reorder Calalnoida and Euphausiacea for plotting
Phy_glom_18s %>%
  ungroup() %>% 
  select(Order) %>% 
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>% 
  unique() %>% 
  as.matrix()->other_orders 

Phy_glom_18s=Phy_glom_18s%>%
  mutate(Order=factor(Order, levels =c("Calanoida","Euphausiacea",other_orders)))


#PC1 Labels
labels_for_PC1=Phy_glom_18s %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))




Phy_glom_18s %>%
  group_by(Order) %>%
  # filter(Order %in% c("Calanoida","Euphausiacea")) %>%
  # filter(prop < 0.001)  %>%
  # bind_rows(.,low_prop_families) %>%
  # filter(prop > 0.01 | Family=="Other") %>% # Create a new data frame with 'other' category
  ggplot(.,aes(x = as.factor(PC1), y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Arcsine-Squareroot\nProportion of Residual Reads", fill = "Order") +
  # ggtitle("Raw Relative Read Abundances By Cycle and Family") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_18s) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)-> majority_p 


majority_p

Phy_glom_18s %>%
  group_by(Order) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  # filter(prop < 0.001)  %>%
  # bind_rows(.,low_prop_families) %>%
  # filter(prop > 0.01 | Family=="Other") %>% # Create a new data frame with 'other' category
  ggplot(.,aes(x = as.factor(PC1), y = asin(sqrt(prop)), fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Proportion of Residual Reads", fill = "Order") +
  # ggtitle("Raw Relative Read Abundances By Cycle and Family") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_18s) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) -> minority_p 
minority_p
stacked_bar_18s=gridExtra::grid.arrange(majority_p,minority_p,ncol=1)


#PNG & PDF Save
ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_18s.png"),
  plot = stacked_bar_18s,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)

ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_18s.pdf"),
  plot = stacked_bar_18s,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)



### HEAT MAPS
# Aggregating data by Order and PC1
agg_data <- Phy_glom_18s %>%
    group_by(Order) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea")))
  group_by(Order, PC1,Sample_ID_short) %>%
  summarise(n_reads = sum(prop)) %>%
  ungroup()

order_abundance <- agg_data %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(n_reads)) %>%
  arrange((total_abundance)) %>%
  pull(Order)

# Reorder the Order factor based on total abundance
agg_data$Order <- factor(agg_data$Order, levels = order_abundance)

# Creating the heatmap using ggplot
ggplot(agg_data, aes(x = as.factor(PC1), y = Order, fill = asin(sqrt((n_reads))))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +  # Adjust color gradient as needed
  theme_minimal() +
  labs(x = "PC1", y = "Order", fill = "n_reads")+
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
