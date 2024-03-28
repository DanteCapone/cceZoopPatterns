#Plot PCR-Bias mitigated data as a function of PC1



# Libraries, Data, Functions, Etc -----------------------------------------


librarian::shelf(tidyverse, googledrive, stringr,here,phyloseq,
                 extrafont, RColorBrewer)



#Add functions for myself
source(("scripts/helpful_functions/treemap_funs_Capone.R"))
source("scripts/helpful_functions/phyloseq_mapping_funs.R")
source("scripts/helpful_functions/general_helper_functions.R")

#Switches
saving=0

#Metadata
##Lat and lon are funky so add back in from metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  select(-X, -Sizefractionmm,max_size) %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)

#Add depth
depths=read.csv(here("data/physical_environmental_data/sample_depths.csv")) %>%
  select(-X) %>%
  mutate(Sample_ID_short=Sample_ID)


#Volume filtered (add to metadata)
volume_filtered=read.csv(here("data/raw_data/biomass/p2107_bt_volume_filtered.csv"))

#Dryweights
dryweights=read.csv("data/raw_data/biomass/dryweights_forzoopmetab.csv")

env_metadata=metadata %>% left_join(.,volume_filtered, by="Sample_ID_short") %>%
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x)

#Predicted proportions
fido_s1=read.csv(here("data/predicted_og/predicted_og_coi_02_26_2024_s1_phy.csv")) %>%
  select(-X)
fido_s2=read.csv(here("data/predicted_og/predicted_og_coi_02_26_2024_s2_phy.csv")) %>%
  select(-X)
fido_s3=read.csv(here("data/predicted_og/predicted_og_coi_02_26_2024_s3_phy.csv")) %>%
  select(-X)

final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 



# Part 1: All Taxa Patterns  ----------------------------------------------------------------


### PART 1: All taxa analysis 
#All merged by cycle
phy_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  group_by(Sample_ID,size) %>%
  summarise(n_reads=sum(n_reads)) %>%
  left_join(.,env_metadata, by="Sample_ID")

#By taxon
phy_taxa_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = str_extract(coord, "(?<=_)[^_]+(?=\\.)"),
         taxa = if_else(is.na(taxa), 'other', taxa))




#Predicted proportions
fido_s1=read.csv(here("data/predicted_og/predicted_og_coi_02_26_2024_s1_phy.csv")) %>%
  select(-X)
fido_s2=read.csv(here("data/predicted_og/predicted_og_coi_02_26_2024_s2_phy.csv")) %>%
  select(-X)
fido_s3=read.csv(here("data/predicted_og/predicted_og_coi_02_26_2024_s3_phy.csv")) %>%
  select(-X)

final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 

#All merged
phy_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  group_by(Sample_ID,size) %>%
  summarise(n_reads=sum(n_reads)) %>%
  left_join(.,env_metadata, by="Sample_ID")

#By species
phy_taxa_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID") %>%
  mutate(taxa=coord)


unique(phy_taxa_pcr$taxa)

## === Plotting === ##
custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  # Add more colors if needed

# Determine the number of unique taxa
num_taxa <- length(unique(phy_taxa_pcr$taxa))
# Choose a suitable palette from RColorBrewer
palette_name <- ifelse(num_taxa <= 8, "Set2", "Set3")  # Example choice, you can adjust as needed
# Generate the color palette
color_palette <- brewer.pal(n = num_taxa, name = palette_name)


labels_for_map=phy_taxa_pcr %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))


#==== Relative Abundances


#Taxa
phy_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR Bias-Mitigated Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Relative Abundance"),
       fill = "Genus") +
  facet_wrap(~size_fraction, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->phy_pcr_props_plot_spp
phy_pcr_props_plot_spp

 
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_genus_props_coi_bar.pdf"), 
    plot = phy_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_genus_props_coi_bar.png"), 
    plot = phy_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

#Relative Abundances


#==== Biomass
#By cycle
phy_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(.,aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(
    title = "PCR Bias Mitigated Calanoid Copepod Biomass",
    x = "Offshore \u2190 PC1 \u2192 Onshore",
    y = expression("Biomass (mg C " ~ m^-2 * ")"),
    fill = "Cycle"
  ) +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=custom_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_biomass_plot


phy_pcr_biomass_plot


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_all_biomass_scaled_all_sites_coi.pdf"), 
    plot = phy_pcr_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }
#




#By taxa
phy_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR Bias Mitigated Proportional Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Genus") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->phy_pcr_biomass_plot
phy_pcr_biomass_plot
# ez_save(phy_pcr_biomass_plot,"plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites")

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_all_genera_biomass_scaled_coi.pdf"), 
    plot = phy_pcr_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_all_genera_biomass_scaled_coi.png"), 
    plot = phy_pcr_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }




# Raw Reads ---------------------------------------------------------------


## === RAW READS: Repeat Analysis with raw/normalized reads === ##
#Read in the data

#Taxa
coi_taxa=read.csv(here("data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv")) %>%
  mutate(Genus = ifelse(Genus == "Genus", Family, Genus)) %>%
  column_to_rownames("Genus") %>% 
  mutate(Hash=X) %>%
  select(-X)



#Predicted proportions
fido_s1_raw=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_taxa_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Genus", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Genus", all = TRUE)%>%
  column_to_rownames("Genus") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_coi_merged_raw





#Metadata
env_metadata_phy=env_metadata %>%
  column_to_rownames("Sample_ID_dot")



#Make phyloseq objects

#coi
# OTU = otu_table(as.matrix(coi_otu), taxa_are_rows = TRUE)
OTU = otu_table(as.matrix(fido_coi_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(env_metadata_phy)
Phy_raw_coi=phyloseq(OTU, TAX, meta)%>%
  phyloseq_transform_to_long(.)

Phy_props_coi=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Genus=asv_code) %>%
  select(-asv_code)

phy_coi=Phy_props_coi 

Phy_norm_coi <- phyloseq_normalize_median(phyloseq(OTU, TAX, meta)) %>%
  #Transform to proportions
  transform_sample_counts(., function(x) x / sum(x) ) %>%
  phyloseq_transform_to_long(.)


#Abundances

#By Taxa
phy_coi %>%
  # Phy_props_coi %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance",
       fill = "Taxa") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->phy_nreads_props_plot_spp
phy_nreads_props_plot_spp

 
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_genus_props_coi_bar.pdf"), 
    plot = phy_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_genus_props_coi_bar.png"), 
    plot = phy_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

#=== Biomass
phy_coi %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 120))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=color_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_nreads_biomass_plot
phy_nreads_biomass_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_all_genera_biomass_scaled_coi.pdf"), 
    plot = phy_nreads_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_all_genera_biomass_scaled_coi.png"), 
    plot = phy_nreads_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}



# Part 2: Calanoid Copepods -----------------------------------------------

#Filter to calanoida
calanoida_taxa_pcr=phy_taxa_pcr %>% mutate(Genus=taxa) %>%
  left_join(.,coi_taxa %>% rownames_to_column("Genus"), by="Genus") %>%
  filter(Order=="Calanoida")


unique(calanoida_taxa_pcr$taxa)

## === Plotting === ##
custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  # Add more colors if needed

# Determine the number of unique taxa
num_taxa <- length(unique(phy_taxa_pcr$taxa))
# Choose a suitable palette from RColorBrewer
palette_name <- ifelse(num_taxa <= 8, "Set1", "Set3")  # Example choice, you can adjust as needed
# Generate the color palette
color_palette <- brewer.pal(n = num_taxa, name = palette_name)


labels_for_map=calanoida_taxa_pcr %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))


#==== Relative Abundances

#Cycle
calanoida_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = cycle)) +
  geom_bar(stat = "identity") +
  labs(title = "PCR-Bias Mitigated Calanoid Copepod Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Relative Abundance"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_pcr_props_plot
calanoida_pcr_props_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_calanoid_props_scaled_cycle.pdf"), 
    plot = calanoida_pcr_props_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

#Taxa
calanoida_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR Bias-Mitigated Calanoid Copepod Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Relative Abundance"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->calanoida_pcr_props_plot_spp
calanoida_pcr_props_plot_spp


 
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_calanoid_props_scaled_spp_coi.pdf"), 
    plot = calanoida_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_calanoid_props_scaled_spp_coi.png"), 
    plot = calanoida_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }




#==== Biomass
### Scatterplot biomass vs. pcr rel abundance
# calanoida_taxa_pcr %>%
#   group_by(Sample_ID) %>%
#   summarize(
#     n_reads = sum(n_reads),
#     biomass_mg_m2 = sum(biomass_mg_m2),
#     across(where(is.numeric), mean)) %>%
#   ggplot(aes(x=(biomass_dry), y=biomass_mg_m2))+
#   geom_point(aes( size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
#   scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
#   scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
#   # geom_smooth(method = "lm", fullrange=TRUE, se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
#   labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "PCR Bias-Mitigated Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
#   ggtitle("Spearman Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance (COI)")+
#   # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.1, label.y = 0.75)+
#   stat_cor(method="spearman", label.x = 0.1, label.y = 0.9)+
#   guides(size = FALSE, fill=FALSE) +
#   # facet_wrap(~offshore_onshore, nrow=3)+
#   # facet_wrap(~size_fraction, nrow=3)+
#   # facet_wrap(~size_fraction, nrow=3)+
#   # geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5, alpha = 0.3) +  # Add 1-to-1 line with modifications
#   theme_classic()+
#   theme(axis.text.x = element_text(hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))


#By cycle
calanoida_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(.,aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(
    title = "PCR Bias Mitigated Calanoid Copepod Biomass",
    x = "Offshore \u2190 PC1 \u2192 Onshore",
    y = expression("Biomass (mg C " ~ m^-2 * ")"),
    fill = "Cycle"
  ) +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=custom_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_pcr_biomass_plot


calanoida_pcr_biomass_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites.pdf"), 
    plot = calanoida_pcr_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }
#


#By taxa
calanoida_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR Bias Mitigated Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Genus") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->calanoida_pcr_biomass_plot
calanoida_pcr_biomass_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites.pdf"), 
    plot = calanoida_pcr_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }





## === RAW READS: Repeat Analysis with raw/normalized reads === ##
#Read in the data

#Make phyloseq objects
calanoida_raw=phy_coi %>% filter(Order=="Calanoida")


#Abundances
calanoida_raw %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance",
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_nreads_props_plot
calanoida_nreads_props_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_reads_calanoid_props_cycle_all_sites.pdf"), 
    plot = calanoida_nreads_props_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }


#By Taxa
calanoida_raw %>%
  # calanoida_props_coi %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance",
       fill = "Taxa") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->calanoida_nreads_props_plot_spp
calanoida_nreads_props_plot_spp

 
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_reads_calanoida_props_spp_all_sites_coi.pdf"), 
    plot = calanoida_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_reads_calanoida_props_spp_all_sites_coi.png"), 
    plot = calanoida_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

#=== Biomass
calanoida_raw %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 120))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=custom_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_nreads_biomass_plot
calanoida_nreads_biomass_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_reads_calanoid_biomass_scaled_cycle_all_sites.pdf"), 
    plot = calanoida_nreads_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}


calanoida_raw %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 120))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values=color_palette)->calanoida_nreads_biomass_plot_fam
calanoida_nreads_biomass_plot_fam

saving=0
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw_reads_calanoid_biomass_Genus_all_sites.pdf"), 
    plot = calanoida_nreads_biomass_plot_fam,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}





## =======================================



# PCR vs Raw --------------------------------------------------------------



## =========Plot anomalies in relative abundances

#Join PCR and RRA dataframes
calanoida_taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction,PC1,cycle) %>%
  summarise(n_reads_pcr=sum(n_reads))->pcr_join

calanoida_raw %>%
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)->pcr_and_raw

#Add difference column
pcr_and_raw_coi=pcr_and_raw %>%
  mutate(difference_pcr_raw=n_reads_pcr-n_reads_raw)

#Plot relative abundance differences
pcr_and_raw_coi %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = difference_pcr_raw, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR-Raw",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance (Corrected-Raw)",
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  # coord_cartesian(ylim=c(-0.5,0.8))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> pcr_raw_plot
pcr_raw_plot

saving=0
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/pcr_rra_relative_abundance_diff.pdf"), 
    plot = pcr_raw_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}


# Zooscan  ----------------------------------------------------------------


#===== Zooscan comparison

#Need to modify string category for joining
size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)


#Read data and look at biomass
zooscan_by_sample=read.csv(here("data/Zooscan/zooscan_by_sample_biomass.csv")) %>%
  select(-X) %>%
  mutate(Sample_ID=sample_id) %>%
  mutate(size_fraction = case_when(
    size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
    TRUE ~ NA_real_)) %>%
  group_by(Sample_ID)

#Add biomass sum, calanoid biomass and proportion of calanoid biomass
zooscan_calanoid=zooscan_by_sample %>%
  filter(object_annotation_category=="Calanoida") 


#Relative Abundances
zooscan_relative=read.csv(here("data/Zooscan/zoop_calanoid_by_sample_relative_abundance.csv"))%>%
  mutate(Sample_ID=sample_id) %>%
  mutate(size_fraction = case_when(
    size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
    TRUE ~ NA_real_)) %>%
  group_by(PC1,size_fraction,Sample_ID) %>%
  summarise(relative_abundance_zoo=sum(relative_abundance)) 

## ==== Biomass & Biomass proportions plots === #
labels_for_map=zooscan_calanoid %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

zooscan_calanoid %>%
  filter(size_fraction!=5) %>%
  ggplot(aes(x = as.factor(PC1), y = biomass_prop_taxa, fill = cycle)) +
  geom_bar(stat = "identity", width=0.8) +
  labs(title = "Calanoid Copepod Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoid_biomass_prop_zooscan
calanoid_biomass_prop_zooscan

# #Save
# ggsave(
#   filename = here("plots/methods_comparison/zooscan_calanoid_biomass_proportions.pdf"), 
#   plot = calanoid_biomass_zooscan,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# 
# #Biomass
# zooscan_calanoid %>%
#   filter(size_fraction!=5) %>%
#   ggplot(aes(x = as.factor(PC1), y = dryweight_C_mg_m2_taxa, fill = cycle)) +
#   geom_bar(stat = "identity", width=0.8) +
#   labs(title = "Calanoid Copepod Zooscan Biomass",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = expression("Biomass (mg C " ~ m^-2 * ")"),
#        fill = "Cycle") +
#   facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_fill_manual(values = custom_palette) +
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoid_biomass_zooscan
# calanoid_biomass_zooscan



# PCR vs Raw vs Zooscan ---------------------------------------------------



#========== COMPARE: Make combined dataframe for comparing all 3 methods
zooscan_calanoid %>%
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_coi, by=c("PC1","size_fraction"))%>%
  mutate(difference_pcr_zoo_biomass=n_reads_pcr*biomass_mg_m2*0.37-dryweight_C_mg_m2_taxa,
         difference_raw_zoo_biomass=n_reads_raw*biomass_mg_m2*0.37-dryweight_C_mg_m2_taxa) %>%
  mutate(difference_pcr_zoo_biomass_prop=n_reads_pcr-biomass_prop_taxa,
         difference_raw_zoo_biomass_prop=n_reads_raw-biomass_prop_taxa)->pcr_raw_zoo_coi


#Biomass Differences
pcr_raw_zoo_coi %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = difference_pcr_zoo_biomass, fill = cycle.x)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR-Zooscan",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance (Corrected-Zooscan)",
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  # coord_cartesian(ylim=c(-0.5,0.8))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> pcr_zoo_plot
pcr_zoo_plot
saving=0
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/pcr_zooscan_relative_abundance_diff.pdf"), 
    plot = pcr_zoo_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

pcr_raw_zoo_coi %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = difference_raw_zoo, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Zooscan-Raw",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance (Zooscan-Raw)",
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  # coord_cartesian(ylim=c(-0.5,0.8))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> raw_zoo_plot
raw_zoo_plot

saving=0
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/rra_zooscan_relative_abundance_diff.pdf"), 
    plot = raw_zoo_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}


grid.arrange(pcr_raw_plot,raw_zoo_plot,pcr_zoo_plot, nrow=3, as.table=TRUE)



#Biomass Prop Differences
pcr_raw_zoo_coi %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = difference_pcr_zoo_biomass, fill = cycle.x)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR-Zooscan",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance (Corrected-Zooscan)",
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  # coord_cartesian(ylim=c(-0.5,0.8))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> pcr_zoo_plot
pcr_zoo_plot



#=============Grouped Barplot site abundances
#Format Long
pcr_raw_zoo_coi_long <- pivot_longer(pcr_raw_zoo_coi, 
                                 cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr), 
                                 names_to = "Method", 
                                 values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05*relative_abundance, "xx",NA))



# Now, mark the rows where relative_abundance is closer to n_reads_pcr
# pcr_raw_zoo_coi_long$is_closer_to_pcr <- ifelse(pcr_raw_zoo_coi_long$Method == "biomass_prop_taxa" & pcr_raw_zoo_coi_long$diff_pcr < pcr_raw_zoo_coi_long$diff_raw, "*", "")
# pcr_raw_zoo_coi_long$is_closer_to_raw <- ifelse(pcr_raw_zoo_coi_long$Method == "biomass_prop_taxa" & pcr_raw_zoo_coi_long$diff_raw > pcr_raw_zoo_coi_long$diff_pcr, "x", "")



#Conduct stats for differences
# pcr_raw_zoo_coi_long %>%
# group_by(PC1) %>%
# summarise(p_value_raw_vs_zoo = t.test(relative_abundance[Method == "n_reads_raw"], relative_abundance[Method == "relative_abundance_zoo"])$p.value,
#           p_value_pcr_vs_zoo = t.test(relative_abundance[Method == "n_reads_pcr"], relative_abundance[Method == "relative_abundance_zoo"])$p.value) %>%
# left_join(.,pcr_raw_zoo_coi_long, by="PC1")->pcr_raw_zoo_coi_long_w_stats


# Calculate ANOVA
anova_result <- aov(relative_abundance ~ Method * Sample_ID.y, data = pcr_raw_zoo_coi_long)
summary(anova_result)

# Perform Tukey's HSD post hoc test
tukey_result <- TukeyHSD(anova_result)

# Print the results
print(tukey_result)

#Plot grouped bar plot
pcr_raw_zoo_coi_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_point(position = position_dodge(width = 0.9), aes(shape = Method), size = 3) +  # Add points
  # geom_line(position = position_dodge(width = 0.9), aes(group = Method, color = Method), size = 2) +  # Add lines with colors
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences COI",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Biomass ",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#70BF41", "#4F86F7", "#F78D4F"),
                    labels = c("Zooscan", "PCR-corrected", "Raw Reads")) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+# Define linetypes for each Method
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F")) -> grouped_bar_all  # Define colors for each Method

grouped_bar_all

ggsave(
  filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_coi.pdf"),
  # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
  plot = grouped_bar_all,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_coi.png"),
  # filename = here("plots/methods_comparison/line_relative_abundance_diff.png"),
  plot = grouped_bar_all,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)



# Correlations ------------------------------------------------------------

## Correlation plot

# Check distribution
pcr_raw_zoo_coi %>%
  ggplot(., aes(x = asin(sqrt(n_reads_pcr)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.1, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

pcr_raw_zoo_coi %>%
  ggplot(., aes(x = asin(sqrt(n_reads_raw)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.1, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

pcr_raw_zoo_coi %>%
  ggplot(., aes(x = asin(sqrt(biomass_prop_taxa)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels



# Correlations between Zoo-PB, PCR-RA and RRA -----------------------------

#ADd clusters to df for plotting groups
clusters=read.csv(here("data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()

pcr_raw_zoo_coi=pcr_raw_zoo_coi %>%
  left_join(.,clusters,by="PC1")


custom_palette <- c("#5BA3D5", "#66CC66", "#FF4C38", "#FF8F66", "#A085D9")




#Zoo vs. pcr
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  # scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", fullrange=TRUE, se = TRUE, formula = y ~ x, aes(group=offshore_onshore, color=offshore_onshore)) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "PCR Bias-Mitigated Relative Abundance\n (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance (COI)")+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.1, label.y = 0.75)+
  stat_cor(method="pearson", label.x = 0.1, label.y = 0.9)+
  guides(size = FALSE, fill=FALSE) +
  # facet_wrap(~offshore_onshore, nrow=3)+
  # facet_wrap(~size_fraction, nrow=3)+
  # facet_wrap(~size_fraction, nrow=3)+
  # geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5, alpha = 0.3) +  # Add 1-to-1 line with modifications
  theme_classic()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +  # Set x limits and ticks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+  # Set y limits and ticks
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr
zoo_vs_pcr
ggsave(
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi_cluster.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi_cycle.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi_size.pdf"),
  filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi.pdf"),
  plot = zoo_vs_pcr,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi_cluster.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi_cycle.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi_size.png"),
  filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_coi.png"),
  plot = zoo_vs_pcr,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)


#### ===== 
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x, fullrange=TRUE) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Read Relative Abundance\n (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion\n and Raw Read Relative Abundance (COI)")+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.1, label.y = 0.9)+
  stat_cor(method="pearson", label.x = 0.1, label.y = 0.9)+
  guides(size = FALSE, fill=FALSE) +
  # facet_wrap(~offshore_onshore, nrow=3)+
  # facet_wrap(~cycle.y, nrow=3)+
  # facet_wrap(~size_fraction, nrow=3)+
  # geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5, alpha = 0.3) +  # Add 1-to-1 line with modifications
  theme_classic()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +  # Set x limits and ticks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+  # Set y limits and ticks  geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5, alpha = 0.5, linetype = "dashed") +  # Add 1-to-1 line with modifications
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw
ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_cluster.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_cycle.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_size.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi.pdf"),
    plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_cluster.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_cycle.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_size.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)



### PCR vs Raw
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=n_reads_raw, y=n_reads_pcr))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between PCR Bias-Mitigated and Raw Read Relative Abundance")+
  stat_cor(method="pearson", label.x = 0.1, label.y = 0.9)+
  guides(size = FALSE, fill=FALSE) +
  # geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5, alpha = 0.3) +  # Add 1-to-1 line with modifications
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw
if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_raw_correlation.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
) }


# ANCOVA ------------------------------------------------------------------



## ANCOVA with groupings
library(rstatix)

# ANCOVA with one continuous covariate (e.g., n_reads_raw)


#Raw vs. Zoo
ancova_result_coi_raw <- pcr_raw_zoo_coi %>%
  filter(!is.na(n_reads_raw))%>%
  mutate(t1 = ifelse(cycle.y == "T1", cycle.y, "other")) %>%
  lm(asin(sqrt(n_reads_raw)) ~ asin(sqrt(biomass_prop_taxa))+offshore_onshore, data = .)
#Can remove offshore_onshore and interactions for cycle
anova(ancova_result_coi_raw)
# Diagnostic plots
par(mfrow=c(2,2)) # Create a 2x2 layout for the plots
plot(ancova_result) # Plot diagnostic plots



#==PCR-RA vs. ZooPB ANCOVA
ancova_result_coi_pcr <- pcr_raw_zoo_coi %>%
  filter(!is.na(n_reads_raw))%>%
  mutate(t1 = ifelse(cycle.y == "T1", cycle.y, "other")) %>%
  lm(asin(sqrt(n_reads_pcr)) ~ asin(sqrt(biomass_prop_taxa))*size_fraction*t1+cycle.y, data = .)
#Can remove offshore_onshore and interactions for cycle
anova(ancova_result_coi_pcr)
# Diagnostic plots
par(mfrow=c(2,2)) # Create a 2x2 layout for the plots
plot(ancova_result_coi_pcr) # Plot diagnostic plots





















### ====== RANDOM ANALYSES


## Correlate biomasses
pcr_raw_zoo_18s%>%
  filter(dryweight_C_mg_sum_sample<20) %>% 
  ggplot(aes(x=(biomass_mg_m2), y=dryweight_C_mg_sum_sample))+
  geom_point(aes( size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  geom_smooth(method = "lm", fullrange=TRUE, se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance (COI)")+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.1, label.y = 0.75)+
  stat_cor(method="pearson", label.x = 0.1, label.y = 0.9)+
  guides(size = FALSE, fill=FALSE) +
  # facet_wrap(~offshore_onshore, nrow=3)+
  # facet_wrap(~size_fraction, nrow=3)+
  # facet_wrap(~size_fraction, nrow=3)+
  # geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5, alpha = 0.3) +  # Add 1-to-1 line with modifications
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))





###SCRAP
#Patterns comparison
## Patterns for only where pcr-bias itigated is better
pcr_raw_zoo_coi_long %>%
  filter(is_closer!="x" | Method=="n_reads_pcr") %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  labs(title = "Methods Differences",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Proportion Reads or Biomass ",
       fill = "Method")+
  coord_cartesian(ylim=c(0,1.1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = c( "#70BF41","#4F86F7", "#F78D4F"),
                    labels = c("Zooscan","PCR-corrected","Raw Reads")) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->grouped_bar_all

grouped_bar_all


## Patterns for only where raw itigated is better
pcr_raw_zoo_coi_long %>%
  filter(is_closer!="*" | Method=="n_reads_pcr") %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  labs(title = "Methods Differences",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Proportion Reads or Biomass ",
       fill = "Method")+
  coord_cartesian(ylim=c(0,1.1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = c( "#70BF41","#4F86F7", "#F78D4F"),
                    labels = c("Zooscan","PCR-corrected","Raw Reads")) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->grouped_bar_all

grouped_bar_all



### Finally correlate all with PC1

pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=PC1, y=biomass_prop_taxa))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  facet_wrap(~size_fraction, nrow=3)+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  geom_smooth(method = "lm", fullrange=TRUE, se = TRUE, color = "black", formula = y ~ x)+  # Add linear regression line
  labs(x = "PC1", y = "PCR Bias-Mitigated Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.75)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

