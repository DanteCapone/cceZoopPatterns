#Plot PCR-Bias mitigated data as a function of PC1
librarian::shelf(tidyverse, googledrive, stringr,here,gridextra,phyloseq,
                 extrafont, RColorBrewer)



#Add functions for myself
source(("scripts/helpful_functions/treemap_funs_Capone.R"))
source("scripts/helpful_functions/phyloseq_mapping_funs.R")
source("scripts/helpful_functions/general_helper_functions.R")

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
dryweights=read.csv("data/raw_data/biomass/dryweights_forzoopmetab.csv") %>%
  mutate(biomass_dry=8/3*biomass_dry)

env_metadata=metadata %>% left_join(.,volume_filtered, by="Sample_ID_short") %>%
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x)


# PCR-Bias Mitigated All --------------------------------------------------

#Predicted proportions
fido_s1=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s1_phy.csv")) %>%
  select(-X)
fido_s2=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s2_phy.csv")) %>%
  select(-X)
fido_s3=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s3_phy.csv")) %>%
  select(-X)

final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 



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
fido_s1=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s1_phy.csv")) %>%
  select(-X)
fido_s2=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s2_phy.csv")) %>%
  select(-X)
fido_s3=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s3_phy.csv")) %>%
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
palette_name <- ifelse(num_taxa <= 8, "Set1", "Set3")  # Example choice, you can adjust as needed
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
  labs(title = "PCR Bias-Mitigated Relative Abundance (18S)",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Relative Abundance"),
       fill = "Family") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
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
    filename = here("plots/methods_comparison/pcr_bias_mitigated/PCR_calanoid_props_scaled_spp_18s.pdf"), 
    plot = phy_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/pcr_bias_mitigated/PCR_calanoid_props_scaled_spp_18s.png"), 
    plot = phy_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }




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
  # coord_cartesian(ylim = c(0, 100))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=custom_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_biomass_plot


phy_pcr_biomass_plot

saving=0
if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/pcr_bias_mitigated/PCR_all_biomass_scaled_all_sites_18s.pdf"), 
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
  labs(title = "PCR Bias Mitigated Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Family") +
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



if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/pcr_bias_mitigated/PCR_all_biomass_scaled_all_sites.pdf"), 
  plot = phy_pcr_biomass_plot,
  width = 8,  # Width in inches
  height = 6  # Height in inches
) }

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/pcr_bias_mitigated/PCR_all_biomass_scaled_all_sites.png"), 
    plot = phy_pcr_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }





# Raw Reads ---------------------------------------------------------------

## === RAW READS: Repeat Analysis with raw/normalized reads === ##
#Read in the data

#Taxa
zhan_taxa=read.csv(here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv")) %>%
column_to_rownames("Family") %>% 
  mutate(Hash=X) %>%
  select(-X, -Species, -Genus)


#Predicted proportions
fido_s1_raw=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Family", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Family", all = TRUE)%>%
  column_to_rownames("Family") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_raw





#Metadata
env_metadata_phy=env_metadata %>%
  column_to_rownames("Sample_ID_dot")



#Make phyloseq objects

#18s
# OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(env_metadata_phy)
Phy_raw_18s=phyloseq(OTU, TAX, meta)%>%
  phyloseq_transform_to_long(.)

Phy_props_18s=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>%
  select(-asv_code)

phy_18s=Phy_props_18s 



#By Taxa
phy_18s %>%
  # Phy_props_18s %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = Family)) +
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
  filename = here("plots/methods_comparison/raw/raw_reads_phy_props_family_all_sites.pdf"), 
  plot = phy_nreads_props_plot_spp,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)}

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw/raw_reads_phy_props_family_all_sites.png"), 
    plot = phy_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

#=== Biomass
phy_18s %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Carbon Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Family") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 120))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values = color_palette)->phy_nreads_biomass_plot
phy_nreads_biomass_plot


if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/raw/raw_reads_biomass_family.pdf"), 
  plot = phy_nreads_biomass_plot,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)}

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw/raw_reads_biomass_family.png"), 
    plot = phy_nreads_biomass_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}



### ============== PART 2: Calanoids ==============

#Filter to calanoida
calanoida_taxa_pcr=phy_taxa_pcr %>% mutate(Family=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Family"), by="Family") %>%
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


labels_for_map=calanoida %>% 
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
    filename = here("plots/methods_comparison/PCR_calanoid_rela_abundance_cycle.pdf"), 
    plot = calanoida_pcr_props_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

#Taxa
calanoida_taxa_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR Bias-Mitigated Calanoid Copepod Relative Abundance (18S)",
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

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/PCR_calanoid_props_scaled_spp_18s.pdf"), 
    plot = calanoida_pcr_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }




#==== Biomass
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
  # coord_cartesian(ylim = c(0, 100))+
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
       fill = "Family") +
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

#Taxa
zhan_taxa=read.csv(here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv")) %>%
  column_to_rownames("Family") %>% 
  mutate(Hash=X) %>%
  select(-X, -Species, -Genus)


#Predicted proportions
fido_s1_raw=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Family", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Family", all = TRUE)%>%
  column_to_rownames("Family") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_raw





#Metadata
env_metadata_phy=env_metadata %>%
  column_to_rownames("Sample_ID_dot")



#Make phyloseq objects
calanoida_raw=phy_18s %>% filter(Order=="Calanoida")


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
  # calanoida_props_18s %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Relative Abundance (18S)",
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

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw/raw_reads_calanoida_props_spp_all_sites_18S.pdf"), 
    plot = calanoida_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}
if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw/raw_reads_calanoida_props_spp_all_sites_18S.png"), 
    plot = calanoida_nreads_props_plot_spp,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}

#=== Biomass


calanoida_raw %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  # coord_cartesian(ylim = c(0, 120))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  scale_fill_manual(values=color_palette)->calanoida_nreads_biomass_plot_fam
calanoida_nreads_biomass_plot_fam

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/raw/raw_reads_calanoid_biomass_Family_all_sites_18s.pdf"), 
    plot = calanoida_nreads_biomass_plot_fam,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}





## =======================================




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
pcr_and_raw_18s=pcr_and_raw %>%
  mutate(difference_pcr_raw=n_reads_pcr-n_reads_raw)

#Plot relative abundance differences
pcr_and_raw_18s %>%
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

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/pcr_rra_relative_abundance_diff.pdf"), 
    plot = pcr_raw_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}


# Zooscan -----------------------------------------------------------------


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
labels_for_map=biomass_map %>% 
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

#Save
ggsave(
  filename = here("plots/methods_comparison/zooscan_calanoid_biomass_proportions.pdf"), 
  plot = calanoid_biomass_zooscan,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)


#Biomass
zooscan_calanoid %>%
  filter(size_fraction!=5) %>%
  ggplot(aes(x = as.factor(PC1), y = dryweight_C_mg_m2_taxa, fill = cycle)) +
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
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoid_biomass_zooscan
calanoid_biomass_zooscan

#Save
ggsave(
  filename = here("plots/methods_comparison/zooscan_calanoid_biomass_proportions.pdf"), 
  plot = calanoid_biomass_zooscan,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)


#========== COMPARE: Make combined dataframe for comparing all 3 methods
zooscan_calanoid %>%
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_18s, by=c("PC1","size_fraction"))%>%
  mutate(difference_pcr_zoo_biomass=n_reads_pcr*biomass_mg_m2*0.37-dryweight_C_mg_m2_taxa,
         difference_raw_zoo_biomass=n_reads_raw*biomass_mg_m2*0.37-dryweight_C_mg_m2_taxa) %>%
  mutate(difference_pcr_zoo_biomass_prop=n_reads_pcr-biomass_prop_taxa,
       difference_raw_zoo_biomass_prop=n_reads_raw-biomass_prop_taxa)->pcr_raw_zoo_18s
  

#Biomass Differences
pcr_raw_zoo_18s %>%
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

if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/pcr_zooscan_relative_abundance_diff.pdf"), 
  plot = pcr_zoo_plot,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)}

pcr_raw_zoo_18s %>%
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


if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/rra_zooscan_relative_abundance_diff.pdf"), 
  plot = raw_zoo_plot,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)}


grid.arrange(pcr_raw_plot,raw_zoo_plot,pcr_zoo_plot, nrow=3, as.table=TRUE)



#Biomass Prop Differences
pcr_raw_zoo_18s %>%
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
pcr_raw_zoo_18s_long <- pivot_longer(pcr_raw_zoo_18s, 
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
  



# ANOVA test
anova_result <- aov(relative_abundance ~ Method * as.factor(PC1), data = pcr_raw_zoo_18s_long)
anova_summary <- summary(anova_result)
print(anova_summary)

# Kruskal-Wallis test (non-parametric alternative)
kruskal_result <- pcr_raw_zoo_18s_long %>%
  group_by(Sample_ID.y) %>%
  do(kruskal_test = kruskal.test(relative_abundance ~ Method, data = .))

# Print results
print(kruskal_result$kruskal_test)

# Perform Tukey's HSD post hoc test
tukey_result <- TukeyHSD(anova_result)

# Print the results
print(tukey_result)

#Plot grouped bar plot
pcr_raw_zoo_18s_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_point(position = position_dodge(width = 0.9), aes(shape = Method), size = 3) +  # Add points
  # geom_line(position = position_dodge(width = 0.9), aes(group = Method, color = Method), size = 2) +  # Add lines with colors
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences",
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
  filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s.pdf"),
  # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s.png"),
  # filename = here("plots/methods_comparison/line_relative_abundance_diff.png"),
  plot = grouped_bar_all,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)
 

## Correlation plot

# Check distribution
pcr_raw_zoo_18s %>%
ggplot(., aes(x = asin(sqrt(n_reads_pcr)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

pcr_raw_zoo_18s %>%
  ggplot(., aes(x = asin(sqrt(n_reads_raw)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

pcr_raw_zoo_18s %>%
  ggplot(., aes(x = asin(sqrt(biomass_prop_taxa)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels



# Correlation
custom_palette <- c("#5BA3D5", "#66CC66", "#FF4C38", "#FF8F66", "#A085D9")

#ADd clusters to df for plotting groups
clusters=read.csv(here("data/physical_environmental_data/pca_clusters.csv")) %>%
                    select(-Sample_ID_dot) %>% unique()
                  
                  
pcr_raw_zoo_18s=pcr_raw_zoo_18s %>%
left_join(.,clusters,by="PC1")
#Check outliers
# Calculate quantiles
quantiles <- quantile(pcr_raw_zoo_18s$biomass_prop_taxa, c(0.25, 0.75))
IQR <- quantiles[2] - quantiles[1]

# Define lower and upper bounds (e.g., using 1.5*IQR)
lower_bound <- quantiles[1] - 1.5 * IQR
upper_bound <- quantiles[2] + 1.5 * IQR
pcr_raw_zoo_18s_outliers_rm=pcr_raw_zoo_18s %>%
  filter(biomass_prop_taxa >= lower_bound & biomass_prop_taxa <= upper_bound)

pcr_raw_zoo_18s_outliers_rm %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "PCR Bias-Mitigated Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Spearman Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "spearman", label.x = 0.1, label.y = 1.3)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr
  zoo_vs_pcr
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_size.pdf"),
    # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_cycle.pdf"),
    # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clust.pdf"),
    # filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_sig_diffs.pdf"), 
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_size.png"),
    # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_cycle.png"),
    # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clust.png"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  

#### =====  
  pcr_raw_zoo_18s_outliers_rm %>%
    filter(!is.na(cycle.y))%>%
    # filter(cycle.y=="1") %>%
    ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_raw))))+
    geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
    scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
    scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
    # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
    labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
    ggtitle("Spearman Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
    stat_cor(method = "spearman", label.x = 0.1, label.y = 1.5)+
    guides(size = FALSE, fill=FALSE) +
    facet_wrap(~size_fraction, nrow=3) +
    theme_classic()+
    theme(axis.text.x = element_text(hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14))->zoo_vs_raw
  zoo_vs_raw
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_size.pdf"),
    # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_cycle.pdf"),
    # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clust.pdf"),
    # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_cycle.pdf"),
    plot = zoo_vs_raw,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )

  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_size.png"),
    # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_cycle.png"),
    # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clust.png"),
    # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clust.png"),
    
    plot = zoo_vs_raw,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  
  
  ### PCR vs Raw
  pcr_raw_zoo_18s_outliers_rm %>%
    filter(!is.na(cycle.y))%>%
    # filter(cycle.y=="1") %>%
    ggplot(.,aes(x=n_reads_pcr, y=n_reads_raw))+
    geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
    scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
    scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
    geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
    labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
    ggtitle("Spearman Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
    stat_cor(method = "spearman", label.x = 0.1, label.y = 0.2)+
    guides(size = FALSE, fill=FALSE) +
    theme_classic()+
    theme(axis.text.x = element_text(hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14))->zoo_vs_raw
  zoo_vs_raw
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation.pdf"),
    plot = zoo_vs_raw,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation.png"),
    plot = zoo_vs_raw,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )


  
  # ANCOVA ------------------------------------------------------------------
  
  
  
  ## ANCOVA with groupings
  library(rstatix)
  #18s raw-RA
  ancova_result_18s_raw <- pcr_raw_zoo_18s_outliers_rm %>%
    filter(!is.na(n_reads_raw))%>%
    mutate(t1 = ifelse(cycle.y == "T1", cycle.y, "other")) %>%
    lm(asin(sqrt(n_reads_raw)) ~ asin(sqrt(biomass_prop_taxa))+size_fraction, data = .)
  
  #Can remove offshore_onshore and interactions for cycle
  anova(ancova_result_18s_raw)
  # Diagnostic plots
  par(mfrow=c(2,2)) # Create a 2x2 layout for the plots
  plot(ancova_result) # Plot diagnostic plots
  
  #18s PCR-RA
  ancova_result_18s_pcr <- pcr_raw_zoo_18s_outliers_rm %>%
    filter(!is.na(n_reads_raw))%>%
    mutate(t1 = ifelse(cycle.y == "T1", cycle.y, "other")) %>%
    lm(asin(sqrt(n_reads_pcr)) ~ asin(sqrt(biomass_prop_taxa))*size_fraction*cycle.y*offshore_onshore, data = .)
  
  #Can remove offshore_onshore and interactions for cycle
  anova(ancova_result_18s_pcr)
  # Diagnostic plots
  par(mfrow=c(2,2)) # Create a 2x2 layout for the plots
  plot(ancova_result) # Plot diagnostic plots




### ========= Scrapped =========== ###
  
  
  # #Patterns comparison
  # ## Patterns for only where pcr-bias itigated is better
  # pcr_raw_zoo_18s_long %>%
  #   filter(is_closer!="x" | Method=="n_reads_pcr") %>%
  #   ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  #   geom_bar(stat = "identity", position = "dodge")+
  #   facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  #   theme_minimal()+
  #   labs(title = "Methods Differences",
  #        x = "Offfshore \u2190 PC1 \u2192 Onshore",
  #        y = "Proportion Reads or Biomass ",
  #        fill = "Method")+
  #   coord_cartesian(ylim=c(0,1.1))+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
  #         axis.text.y = element_text(size = 12),
  #         axis.title = element_text(size = 14),
  #         strip.text = element_text(size = 14))+
  #   scale_fill_manual(values = c( "#70BF41","#4F86F7", "#F78D4F"),
  #                     labels = c("Zooscan","PCR-corrected","Raw Reads")) +
  #   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->grouped_bar_all
  # 
  # grouped_bar_all
  # 
  # 
  # ## Patterns for only where raw itigated is better
  # pcr_raw_zoo_18s_long %>%
  #   filter(is_closer!="*" | Method=="n_reads_pcr") %>%
  #   ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  #   geom_bar(stat = "identity", position = "dodge")+
  #   facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  #   theme_minimal()+
  #   labs(title = "Methods Differences",
  #        x = "Offfshore \u2190 PC1 \u2192 Onshore",
  #        y = "Proportion Reads or Biomass ",
  #        fill = "Method")+
  #   coord_cartesian(ylim=c(0,1.1))+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
  #         axis.text.y = element_text(size = 12),
  #         axis.title = element_text(size = 14),
  #         strip.text = element_text(size = 14))+
  #   scale_fill_manual(values = c( "#70BF41","#4F86F7", "#F78D4F"),
  #                     labels = c("Zooscan","PCR-corrected","Raw Reads")) +
  #   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->grouped_bar_all
  # 
  # grouped_bar_all
  
  
  
  
  
# 
# 
# ## === Plotting === ##
# custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  # Add more colors if needed
# 
# labels_for_map=calanoida %>% 
#   ungroup()%>%
#   select(Sample_ID_short,PC1) %>%
#   unique(.) %>%
#   arrange((PC1))
# 
# #Biomass
# #By cycle
# calanoida_taxa_pcr %>%
#   # filter(biomass_mg_m2 >0) %>%
#   ggplot(.,aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(
#     title = "PCR Bias Mitigated Calanoid Copepod Biomass",
#     x = "Offshore \u2190 PC1 \u2192 Onshore",
#     y = expression("Biomass (mg C " ~ m^-2 * ")"),
#     fill = "Cycle"
#   ) +
#   facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   coord_cartesian(ylim = c(0, 100))+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_fill_manual(values=custom_palette)+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_biomass_plot
# 
# 
# phy_pcr_biomass_plot
# ggsave(
#   filename = here("plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites.pdf"), 
#   plot = phy_pcr_biomass_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# #
# 
# 
# #By taxa
# calanoida_taxa_pcr %>%
#   # filter(biomass_mg_m2 >0) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = taxa)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "PCR Bias Mitigated Calanoid Copepod Biomass",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = expression("Biomass (mg C " ~ m^-2 * ")"),
#        fill = "Cycle") +
#   facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
#   coord_cartesian(ylim = c(0, 100))+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_biomass_plot
# phy_pcr_biomass_plot
# # ez_save(phy_pcr_biomass_plot,"plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites")
# ggsave(
#   filename = here("plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites.pdf"), 
#   plot = phy_pcr_biomass_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# 
# #By Species
# phy_spp %>%
#   # filter(biomass_mg_m2 >0) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads*(biomass_mg_m2)*0.37, fill = taxa)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "PCR Bias Mitigated Calanoid Copepod Biomass",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = expression("Biomass (mg C " ~ m^-2 * ")"),
#        fill = "Cycle") +
#   facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_biomass_plot_spp
# phy_pcr_biomass_plot_spp
# ggsave(
#   filename = here("plots/methods_comparison/PCR_calanoid_biomass_scaled_spp.pdf"), 
#   plot = phy_pcr_biomass_plot_spp,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# #Relative Abundances
# 
# #Cycle
# calanoida_taxa_pcr %>%
#   # filter(biomass_mg_m2 >0) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads, fill = cycle)) +
#   geom_bar(stat = "identity") +
#   labs(title = "PCR-Bias Mitigated Calanoid Copepod Relative Abundance",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = expression("Relative Abundance"),
#        fill = "Cycle") +
#   facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_props_plot
# phy_pcr_props_plot
# ggsave(
#   filename = here("plots/methods_comparison/PCR_calanoid_props_scaled_cycle.pdf"), 
#   plot = phy_pcr_props_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# #Taxa
# calanoida_taxa_pcr %>%
#   # filter(biomass_mg_m2 >0) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads, fill = taxa)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "PCR Bias-Mitigated Calanoid Copepod Relative Abundance",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = expression("Relative Abundance"),
#        fill = "Cycle") +
#   facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_pcr_props_plot_spp
# phy_pcr_props_plot_spp
# 
# 
# ggsave(
#   filename = here("plots/methods_comparison/PCR_calanoid_props_scaled_spp_all_sites.pdf"), 
#   plot = phy_pcr_props_plot_spp,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# ## === RAW READS: Repeat Analysis with raw/normalized reads === ##
# #Read in the data
# #Predicted proportions
# fido_s1_raw=read.csv(here("data/fido/fido_18s_s1_ecdf_hash.csv")) %>% 
#   select(-starts_with("X")) %>% 
#   pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
#   mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
#   group_by(Sample_ID_short, Hash) %>%
#   summarise(n_reads = sum(n_reads)) %>%
#   filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
#   pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)
# 
# 
# fido_s2_raw=read.csv(here("data/fido/fido_18s_s2_ecdf_hash.csv")) %>% 
#   select(-starts_with("X")) %>% 
#   pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
#   mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
#   group_by(Sample_ID_short, Hash) %>%
#   summarise(n_reads = sum(n_reads)) %>%
#   filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
#   pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)
# 
# 
# fido_s3_raw=read.csv(here("data/fido/fido_18s_s3_ecdf_hash.csv")) %>% 
#   select(-starts_with("X")) %>% 
#   pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
#   mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
#   group_by(Sample_ID_short, Hash) %>%
#   summarise(n_reads = sum(n_reads)) %>%
#   filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
#   pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)
# 
# 
# merge(fido_s1_raw, fido_s2_raw, by = "Hash", all = TRUE) %>%
#   merge(.,fido_s3_raw, by = "Hash", all = TRUE)%>%
#   column_to_rownames("Hash") %>%
#   mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_raw
# 
# # a=fido_18s_merged_raw %>%
# #   rownames_to_column("Hash") %>%
# #   left_join(zhan_taxa)
# #18s
# # zhan_otu=read.csv(here("data/raw_reads/otu_18s_prefilt_all.csv"), header=TRUE, row.names = 1) %>%
# #   select(where(~ !any(is.na(.)))) %>%
# #   #Convert to proportion
# #   mutate(across(everything(), ~ . / sum(.)))
# 
# 
# zhan_taxa=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_tax.csv")) %>%
#   column_to_rownames("Hash")
# 
# #Metadata
# env_metadata_phy=env_metadata %>%
#   column_to_rownames("Sample_ID_dot")
# 
# 
# 
# #Make phyloseq objects
# 
# #18s
# # OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
# OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
# TAX = tax_table(as.matrix(zhan_taxa))
# meta=sample_data(env_metadata_phy)
# Phy_raw_18s=phyloseq(OTU, TAX, meta)%>%
#   phyloseq_transform_to_long(.)
# 
# Phy_props_18s=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
#   phyloseq_transform_to_long(.)%>%
#   mutate(Order = if_else(is.na(Order), Class, Order)) %>%
#   mutate(Order = if_else(Order=="", Class, Order)) %>%
#   
#   
#   mutate(Family = if_else(is.na(Family), Order, Family)) %>%
#   mutate(Family = if_else(Family=="", Order, Family)) %>%
#   
#   mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
#   mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
#   
#   mutate(Species = if_else(is.na(Species), Genus, Species))%>%
#   mutate(Species = if_else(Species== "", Genus, Species))
# 
# phy_18s=Phy_props_18s %>%
#   filter(Order=="Calanoida")
# 
# Phy_norm_18s <- phyloseq_normalize_median(phyloseq(OTU, TAX, meta)) %>%
#   #Transform to proportions
#   transform_sample_counts(., function(x) x / sum(x) ) %>%
#   phyloseq_transform_to_long(.)
# 
# 
# #Abundances
# phy_18s %>%
#   # filter(Species %in% phy_spp$taxa) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads, fill = cycle)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "Raw Reads Calanoid Copepod Relative Abundance",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = "Relative Abundance",
#        fill = "Cycle") +
#   facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_nreads_props_plot
# phy_nreads_props_plot
# 
# ggsave(
#   filename = here("plots/methods_comparison/raw_reads_calanoid_props_cycle_all_sites.pdf"), 
#   plot = phy_nreads_props_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# 
# #By Taxa
# phy_18s %>%
#   # Phy_props_18s %>%
#   # filter(Species %in% phy_spp$taxa) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads, fill = Species)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "Raw Reads Calanoid Copepod Relative Abundance",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = "Relative Abundance",
#        fill = "Taxa") +
#   facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_nreads_props_plot_spp
# phy_nreads_props_plot_spp
# 
# ggsave(
#   filename = here("plots/methods_comparison/raw_reads_calanoid_props_spp_all_sites.pdf"), 
#   plot = phy_nreads_props_plot_spp,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# #=== Biomass
# phy_18s %>%
#   # filter(Species %in% phy_spp$taxa) %>%
#   ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "Raw Reads Calanoid Copepod Biomass",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = expression("Biomass (mg C " ~ m^-2 * ")"),
#        fill = "Cycle") +
#   facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   coord_cartesian(ylim = c(0, 120))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   scale_fill_manual(values=custom_palette)+
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->phy_nreads_biomass_plot
# phy_nreads_biomass_plot
# 
# 
# ggsave(
#   filename = here("plots/methods_comparison/raw_reads_calanoid_biomass_scaled_cycle_all_sites.pdf"), 
#   plot = phy_nreads_biomass_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## =========PART 3: Anomalies anomalies in relative abundances
# 
# #Join PCR and RRA dataframes
# phy_taxa_pcr %>%
#   mutate(Sample_ID=Sample_ID_dot) %>%
#   group_by(Sample_ID,size_fraction,PC1,cycle) %>%
#   summarise(n_reads_pcr=sum(n_reads))->pcr_join
# 
# phy_18s %>%
#   group_by(Sample_ID,size_fraction) %>%
#   summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
#   left_join(pcr_join, by="Sample_ID") %>%
#   select(-size_fraction.x) %>%
#   mutate(size_fraction=size_fraction.y)->pcr_and_raw
# 
# #Add difference column
# pcr_and_raw=pcr_and_raw %>%
#   mutate(difference_pcr_raw=n_reads_pcr-n_reads_raw)
# 
# #Plot relative abundance differences
# pcr_and_raw %>%
#   # filter(Species %in% phy_spp$taxa) %>%
#   ggplot(aes(x = as.factor(PC1), y = difference_pcr_raw, fill = cycle)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "PCR-Raw",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = "Relative Abundance (Corrected-Raw)",
#        fill = "Cycle") +
#   facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   # coord_cartesian(ylim=c(-0.5,0.8))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> pcr_raw_plot
# 
# ggsave(
#   filename = here("plots/methods_comparison/pcr_rra_relative_abundance_diff.pdf"), 
#   plot = pcr_raw_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# 
# #===== Zooscan comparison
# size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)
# 
# zooscan_by_sample=read.csv(here("data/Zooscan/zooscan_by_sample_biomass.csv")) %>%
#   select(-X) %>%
#   mutate(Sample_ID=sample_id) %>%
#   mutate(size_fraction = case_when(
#     size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
#     TRUE ~ NA_real_)) %>%
#   group_by(Sample_ID) %>%
#   mutate(total_biomass = sum(dryweight_C_mg_sum, na.rm = TRUE))
# 
# #Add biomass sum
# 
# zooscan_calanoid=read.csv(here("data/Zooscan/zoop_calanoid_by_sample.csv")) %>%
#   select(-X) %>%
#   mutate(Sample_ID=sample_id) %>%
#   mutate(size_fraction = case_when(
#     size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
#     TRUE ~ NA_real_)) %>%
#   group_by(PC1,size_fraction,Sample_ID) %>%
#   summarise(dryweight_raw=sum(dryweight_C_mg_m2))
# 
# zooscan_relative=read.csv(here("data/Zooscan/zoop_calanoid_by_sample_relative_abundance.csv"))%>%
#   mutate(Sample_ID=sample_id) %>%
#   mutate(size_fraction = case_when(
#     size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
#     TRUE ~ NA_real_)) %>%
#   group_by(PC1,size_fraction,Sample_ID) %>%
#   summarise(relative_abundance_zoo=sum(relative_abundance)) 
# 
# zooscan_relative %>%
#   filter(size_fraction != 5) %>%
#   left_join(pcr_and_raw, by=c("PC1","size_fraction"))%>%
#   mutate(difference_pcr_zoo=n_reads_pcr-relative_abundance_zoo)%>%
#   mutate(difference_raw_zoo=n_reads_raw-relative_abundance_zoo)->pcr_raw_zoo
# 
# 
# 
# pcr_raw_zoo %>%
#   # filter(Species %in% phy_spp$taxa) %>%
#   ggplot(aes(x = as.factor(PC1), y = difference_pcr_zoo, fill = cycle)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "PCR-Zooscan",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = "Relative Abundance (Corrected-Zooscan)",
#        fill = "Cycle") +
#   facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   # coord_cartesian(ylim=c(-0.5,0.8))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> pcr_zoo_plot
# pcr_zoo_plot
# ggsave(
#   filename = here("plots/methods_comparison/pcr_zooscan_relative_abundance_diff.pdf"), 
#   plot = pcr_zoo_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# pcr_raw_zoo %>%
#   # filter(Species %in% phy_spp$taxa) %>%
#   ggplot(aes(x = as.factor(PC1), y = difference_raw_zoo, fill = cycle)) +
#   geom_bar(stat = "identity", position = "stack", width=0.8) +
#   labs(title = "Zooscan-Raw",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = "Relative Abundance (Zooscan-Raw)",
#        fill = "Cycle") +
#   facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   # coord_cartesian(ylim=c(-0.5,0.8))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> raw_zoo_plot
# raw_zoo_plot
# 
# ggsave(
#   filename = here("plots/methods_comparison/rra_zooscan_relative_abundance_diff.pdf"), 
#   plot = raw_zoo_plot,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
# 
# 
# grid.arrange(pcr_raw_plot,raw_zoo_plot,pcr_zoo_plot, nrow=3, as.table=TRUE)
# 
# 
# 
# #Grouped site abundances
# 
# #Format Long
# pcr_raw_zoo_long <- pivot_longer(pcr_raw_zoo, 
#                                  cols = c(relative_abundance_zoo, n_reads_raw, n_reads_pcr), 
#                                  names_to = "Method", 
#                                  values_to = "relative_abundance") %>%
#   filter(!is.na(Sample_ID.y)) %>%
#   mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
#          diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))
# 
# 
# # Now, mark the rows where relative_abundance is closer to n_reads_pcr
# pcr_raw_zoo_long$is_closer_to_pcr <- ifelse(pcr_raw_zoo_long$Method == "relative_abundance_zoo" & pcr_raw_zoo_long$diff_pcr < pcr_raw_zoo_long$diff_raw, "*", "")
# pcr_raw_zoo_long$is_closer_to_raw <- ifelse(pcr_raw_zoo_long$Method == "relative_abundance_zoo" & pcr_raw_zoo_long$diff_raw < pcr_raw_zoo_long$diff_pcr, "x", "")
# 
# 
# 
# #Conduct stats for differences
# # pcr_raw_zoo_long %>%
# # group_by(PC1) %>%
# # summarise(p_value_raw_vs_zoo = t.test(relative_abundance[Method == "n_reads_raw"], relative_abundance[Method == "relative_abundance_zoo"])$p.value,
# #           p_value_pcr_vs_zoo = t.test(relative_abundance[Method == "n_reads_pcr"], relative_abundance[Method == "relative_abundance_zoo"])$p.value) %>%
# # left_join(.,pcr_raw_zoo_long, by="PC1")->pcr_raw_zoo_long_w_stats
# 
# 
# # Calculate ANOVA
# anova_result <- aov(relative_abundance ~ Method * Sample_ID.y, data = pcr_raw_zoo_long)
# summary(anova_result)
# 
# # Perform Tukey's HSD post hoc test
# tukey_result <- TukeyHSD(anova_result)
# 
# # Print the results
# print(tukey_result)
# 
# #Plot grouped bar plot
# ggplot(pcr_raw_zoo_long, aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
#   geom_bar(stat = "identity", position = "dodge")+
#   facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
#   theme_minimal()+
#   labs(title = "Methods Differences",
#        x = "Offfshore \u2190 PC1 \u2192 Onshore",
#        y = "Relative Abundance",
#        fill = "Method")+
#   coord_cartesian(ylim=c(0,1.25))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 14))+
#   geom_text(aes(label = pcr_raw_zoo_long$is_closer_to_pcr), position = position_dodge(width = 0.9), vjust = -0.1, size = 7) +
#   geom_text(aes(label = pcr_raw_zoo_long$is_closer_to_raw), position = position_dodge(width = 0.9), vjust = -0.15, size = 5) +
#   scale_fill_manual(values = c( "#70BF41","#4F86F7", "#F78D4F"), 
#                     labels = c("PCR-corrected", "Raw reads", "Zooscan")) +
#   scale_x_discrete(labels = labels_for_map$Sample_ID_short)->grouped_bar_all
# 
# grouped_bar_all
# 
# ggsave(
#   filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff.pdf"), 
#   plot = grouped_bar_all,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
