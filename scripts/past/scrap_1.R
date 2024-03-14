#Plot PCR-Bias mitigated data as a function of PC1
librarian::shelf(tidyverse, googledrive, stringr,here,gridextra,phyloseq,
                 extrafont)

#Add functions for myself
source(("scripts/helpful_functions/treemap_funs_Capone.R"))
source("scripts/helpful_functions/phyloseq_mapping_funs.R")
source("scripts/helpful_functions/general_helper_functions.R")


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
volume_filtered=read.csv(here("data/biomass/p2107_bt_volume_filtered.csv"))

#Dryweights
dryweights=read.csv("data/biomass/dryweights_forzoopmetab.csv")

env_metadata=metadata %>% left_join(.,volume_filtered, by="Sample_ID_short") %>%
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x)

#Predicted proportions
fido_s1=read.csv(here("data/predicted_og/predicted_og_coi_02_07_2024_s1.csv")) %>%
  select(-X)
fido_s2=read.csv(here("data/predicted_og/predicted_og_coi_02_07_2024_s2.csv")) %>%
  select(-X)
fido_s3=read.csv(here("data/predicted_og/predicted_og_coi_02_07_2024_s3.csv")) %>%
  select(-X)

final_data_all_sizes=read.csv(here("data/predicted_og/predicted_og_coi_02_07_2024_all.csv")) %>%
  select(-X)%>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 

#All merged
calanoida_pcr= final_data_all_sizes %>%
  filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  group_by(Sample_ID,size) %>%
  summarise(n_reads=sum(n_reads)) %>%
  left_join(.,env_metadata, by="Sample_ID")

#By species
calanoida_spp_pcr= final_data_all_sizes %>%
  filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = sub("^.*?\\.(.*)\\..*$", "\\1", coord),
         taxa = ifelse(taxa == "", "Calanoida", taxa)) 


## === Plotting === ##
custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  # Add more colors if needed

labels_for_map=calanoida %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

#Biomass
#By cycle
calanoida_spp_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(
    title = "PCR Bias Mitigated Calanoid Copepod Biomass",
    x = "Offshore \u2190 PC1 \u2192 Onshore",
    y = expression("Biomass (mg C " ~ m^-2 * ")"),
    fill = "Cycle"
  ) +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  coord_cartesian(ylim = c(0, 120))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=custom_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_pcr_biomass_plot


calanoida_pcr_biomass_plot
ez_save(calanoida_pcr_biomass_plot,"plots/methods_comparison/","PCR_calanoid_biomass_scaled_all_sites")
#


#By taxa
calanoida_spp_pcr %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.7) +
  labs(title = "PCR Bias Mitigated Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_pcr_biomass_plot
calanoida_pcr_biomass_plot
ez_save(calanoida_pcr_biomass_plot,"plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites")
# ez_save(calanoida_pcr_biomass_plot,"plots/methods_comparison/PCR_calanoid_biomass_scaled_all_sites")



#By Species
calanoida_spp %>%
  # filter(biomass_mg_m2 >0) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*(biomass_mg_m2)*0.37, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR Bias Mitigated Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_pcr_biomass_plot_spp
calanoida_pcr_biomass_plot_spp


ez_save(calanoida_pcr_biomass_plot_spp,"plots/methods_comparison","PCR_calanoid_biomass_scaled_spp_all_sites")
# ez_save(calanoida_pcr_biomass_plot,"plots/methods_comparison/PCR_calanoid_biomass_scaled_spp")


#Relative Abundances

#Cycle
calanoida_spp_pcr %>%
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
ez_save(calanoida_pcr_props_plot,"plots/methods_comparison","PCR_calanoid_props_scaled_cycle")


#Taxa
calanoida_spp_pcr %>%
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
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_pcr_props_plot_spp
calanoida_pcr_props_plot_spp


ez_save(calanoida_pcr_props_plot_spp,"plots/methods_comparison","PCR_calanoid_props_scaled_spp_all_sites")


## ============== Raw Reads: Repeat Analysis with raw/normalized reads ================= ##
#Read in the data
fido_s1_raw_coi=read.csv(here("data/fido/fido_coi_s1_ecdf.csv"), row.names = 1) %>% 
  select(-starts_with("X")) %>% 
  rownames_to_column("Hash")%>%
  pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Hash) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0) 


fido_s2_raw_coi=read.csv(here("data/fido/fido_coi_s2_ecdf.csv"), row.names = 1) %>% 
  select(-starts_with("X")) %>% 
  rownames_to_column("Hash")%>%
  pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Hash) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw_coi=read.csv(here("data/fido/fido_coi_s3_ecdf.csv"), row.names = 1) %>% 
  select(-starts_with("X")) %>% 
  rownames_to_column("Hash")%>%
  pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Hash) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw_coi, fido_s2_raw_coi, by = "Hash", all = TRUE) %>%
  merge(.,fido_s3_raw_coi, by = "Hash", all = TRUE)%>%
  column_to_rownames("Hash") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_coi_merged_raw


coi_taxa=read.csv(here("data/phyloseq_bio_data/coi/metazooprunedcoi_tax.csv ")) 

# Create a new row with "other" in all columns
other <- data.frame(lapply(coi_taxa, function(x) "other"))

# Bind the new row to the original dataframe
coi_taxa %>% bind_rows(., other) %>%
  column_to_rownames("Hash")-> coi_taxa

#Metadata
env_metadata_phy=metadata %>% left_join(.,volume_filtered, by="Sample_ID_short") %>%
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x)%>%
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
  phyloseq_transform_to_long(.)%>%
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species))

calanoida_coi=Phy_props_coi %>%
  filter(Order=="Calanoida")

Phy_norm_coi <- phyloseq_normalize_median(phyloseq(OTU, TAX, meta)) %>%
  #Transform to proportions
  transform_sample_counts(., function(x) x / sum(x) ) %>%
  phyloseq_transform_to_long(.)


#Abundances
calanoida_coi %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Relative Abundance",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Biomass (mg/m^2)",
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_nreads_props_plot
calanoida_nreads_props_plot

ez_save(calanoida_nreads_props_plot,"plots/methods_comparison/raw_reads_calanoid_props_cycle_all_sites_coi")

#By Taxa
# calanoida_coi %>%
Phy_props_coi %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads, fill = Species)) +
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
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_nreads_props_plot_spp
calanoida_nreads_props_plot_spp

ez_save(calanoida_nreads_props_plot_spp,"plots/methods_comparison/raw_reads_calanoid_props_spp_all_sites_coi")


#=== Biomass
calanoida_coi %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = n_reads*biomass_mg_m2*0.37, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "Raw Reads Calanoid Copepod Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~max_size, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm")))) +
  theme_minimal()+
  coord_cartesian(ylim = c(0, 120))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values=custom_palette)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoida_nreads_biomass_plot
calanoida_nreads_biomass_plot


ez_save(calanoida_nreads_biomass_plot,"plots/methods_comparison/raw_reads_calanoid_biomass_scaled_cycle_all_sites_coi")

## Plot anomalies in relative abundances

#Join PCR and RRA dataframes
calanoida_spp_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction,PC1,cycle) %>%
  summarise(n_reads_pcr=sum(n_reads))->pcr_join

calanoida_coi %>%
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw=sum(n_reads)) %>%  
  left_join(pcr_join, by="Sample_ID")->pcr_and_raw

#Add difference column
pcr_and_raw=pcr_and_raw %>%
  mutate(difference=n_reads_pcr-n_reads_raw)

#Plot relative abundance differences
pcr_and_raw %>%
  # filter(Species %in% calanoida_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = difference, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR-Raw",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Biomass (mg/m^2)",
       fill = "Cycle") +
  facet_wrap(~size_fraction.x, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)
