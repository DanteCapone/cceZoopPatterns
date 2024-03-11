#Zooscan analysis final
librarian::shelf(tidyverse, googledrive, stringr,here,gridextra)
here()
source(here("scripts/Zooscan_analysis/zooscan_functions.R"))



#Load in the data
# List all .tsv files in the folder
tsv_files <- list.files(here("data/Zooscan/"), pattern = "\\.tsv$", full.names = TRUE)

# Find the most recently added file
latest_ecotaxa <- tsv_files[which.max(file.info(tsv_files)$ctime)]
latest_ecotaxa

# Read the .tsv file into a data frame
zooscan_exp <- read.table(latest_ecotaxa, header=TRUE, sep="\t", encoding="latin1")


#Basic plots to show the relative abundance and biomass of zooscan results from each station
zooscan_processed=readEcotaxa(zooscan_exp)%>%
  #Cleaning up fomratting issues
  mutate(sample_id = str_replace_all(sample_id, "-", "_"),
         sample_id = ifelse(sample_id == "c2_t1_h36", "ct2_t1_h36", sample_id),
         sample_id = ifelse(sample_id == "ct1_t8_h10", "c1_t8_h10", sample_id),
         sample_id = ifelse(sample_id == "ct2_t9_h19", "c2_t9_h19", sample_id),
         sample_id = ifelse(sample_id == "c3_bt6_h25", "c3_t6_h25", sample_id)) 


#Calculate C-biomass
 zooscan_processed %>%
   transform_by_taxa_group(.,"esd") %>%
   mutate(dryweight_C_mg=dryweight_C_ug/1000) %>%
   #Add log biomass, and biomass/m2
   mutate(log10_dryweight_C_mg=log10(dryweight_C_mg))%>%
   #Need to fix hyperiids
   filter(object_annotation_category!="Hyperiidea")->zooscan_biomass 

#Biomass histogram
 zooscan_biomass %>% 
   filter(object_annotation_category %in% c("Calanoida","Copepoda<Maxillopoda","Oithonidae","Harpacticoida", "Poecilostomatoida"))%>%
   filter(object_annotation_category=="Calanoida") %>%
 ggplot(., aes(x = log10((dryweight_C_ug)))) +   # Set the data and the variable to plot
   geom_histogram(binwidth = 0.2, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
   labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels


#Compute relative abudances 
relative_abundances=zooscan_processed %>%
  group_by(sample_id,size_fraction,object_annotation_category) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count),
         relative_abundance = count / total)



##Lat and lon are funky so add back in from metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
# Standardize key columns (some samples are incorerctly notated to compare with Zooscan)
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
  select(-X,-Sample_ID_dot, -Sizefractionmm,-max_size) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)

# Merge dataframes: add metadata to relative abudance data
relative_abundances_map <- relative_abundances %>%
  left_join(metadata, by="sample_id") %>%
  distinct(.)

biomass_map = zooscan_biomass %>%
  left_join(metadata, by="sample_id") %>%
  distinct(.) %>%
  #Make size class a factor
  mutate(size_fraction=factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2", ">2")))


#Filter to calalnoids, merge by sample and size and then compute biomass/m2 and log biomass/m2
#General Zooscan dataframe for plotting
zooscan_by_sample = biomass_map %>%
  group_by(size_fraction, sample_id) %>%
  mutate(dryweight_C_mg_sum_sample = sum(dryweight_C_mg, na.rm = TRUE),
            dryweight_C_ug_sum_sample = sum(dryweight_C_ug, na.rm = TRUE)) %>%
  group_by(size_fraction, sample_id,object_annotation_category) %>%
  summarise(dryweight_C_mg_sum_taxa = sum(dryweight_C_mg, na.rm = TRUE),
            dryweight_C_ug_sum_taxa = sum(dryweight_C_ug, na.rm = TRUE),
            dryweight_C_mg_sum_sample=mean(dryweight_C_mg_sum_sample),
            dryweight_C_ug_sum_sample=mean(dryweight_C_ug_sum_sample),
            sample_conc=mean(sample_conc))%>%
  left_join(metadata, by="sample_id") %>%
  mutate(log10_dryweight_C_ug_m2_taxa=log(dryweight_C_ug_sum_taxa),
         dryweight_C_mg_m2_taxa=dryweight_C_mg_sum_taxa*sample_conc,
         log10_dryweight_C_ug_m2_sample=log(dryweight_C_ug_sum_sample),
         dryweight_C_mg_m2_sample=dryweight_C_mg_sum_sample*sample_conc) %>%
  mutate(biomass_prop_taxa=dryweight_C_mg_m2_taxa/dryweight_C_mg_m2_sample)


#CALANOID DATA FRAME FOR PLOTTING
zoop_calanoid_by_sample = zooscan_by_sample %>%
  filter(object_annotation_category=="Calanoida")

#Calanoid dataframe for relative abundances
zoop_calanoid_by_sample_relative_abundance = relative_abundances_map %>%
  filter(object_annotation_category=="Calanoida")

#Save
write.csv(zoop_calanoid_by_sample,here("data/Zooscan/zoop_calanoid_by_sample_biomass.csv"))
write.csv(zoop_calanoid_by_sample_relative_abundance,here("data/Zooscan/zoop_calanoid_by_sample_relative_abundance.csv"))
write.csv(zooscan_by_sample,here("data/Zooscan/zooscan_by_sample_biomass.csv"))


####PLOTTING

#Color pallete
custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  # Add more colors if needed


#==Proportions==#
relative_abundances_map %>%
  filter(object_annotation_category=="Calanoida") %>%
  filter(size_fraction!=">2") %>%
  ggplot(aes(x = as.factor(PC1), y = relative_abundance, fill = cycle)) +
  geom_bar(stat = "identity") +
  labs(title = "Calanoid Copepod Zooscan Relative Abundances",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Relative Abundance"),
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
  theme_minimal()+
  ylim(0, 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoid_props_zooscan
calanoid_props_zooscan

#Save
ez_save(calanoid_props_zooscan,"plots/methods_comparison/zooscan_calanoid_relative_abundances.jpeg")
# ez_save(calanoid_biomass,"plots/Zooscan/zooscan_calanoid_biomass.jpeg")




## ==== biomass === #
labels_for_map=biomass_map %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

zoop_calanoid_by_sample %>%
# biomass_map %>%
#   filter(object_annotation_category=="Calanoida") %>%
  # filter(Sample_ID_short %in% unique(dryweights$Sample_ID_short))
  filter(size_fraction!=">2") %>%
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
# ez_save(calanoid_biomass_zooscan,"plots/methods_comparison/zooscan_calanoid_biomass_scaled.jpeg")
# ez_save(calanoid_biomass,"plots/Zooscan/zooscan_calanoid_biomass.jpeg")
ggsave(
  filename = here("plots/methods_comparison/zooscan_calanoid_biomass_scaled.pdf"), 
  plot = calanoid_biomass_zooscan,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)


#Continuous PC1
biomass_map %>%
  filter(object_annotation_category=="Calanoida") %>%
  ggplot(aes(x = PC1, y = dryweight_C, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack", width=0.2) +
  labs(title = "Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Biomass (ug C)",
       fill = "Taxa Group") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c(">2 mm","0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)









#With continuous PC1
###Zooscan relative abundances vs. pc1
relative_abundances_map %>%
  filter(relative_abundance > 0.05, object_annotation_category != "multiple organisms") %>%
  ggplot(aes(x = PC1, y = relative_abundance, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack", width=0.2) +
  geom_line(aes(x = as.numeric(PC1), y = relative_abundance, color=object_annotation_category), position = position_jitter(width = 0.1, height = 0)) +
  labs(title = "Stacked Bar Plot",
       x = "PC1",
       y = "Relative Abundance",
       fill = "Taxa Group") +
  facet_wrap(~size_fraction, nrow = 4) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10)) +
  scale_x_continuous(labels = function(x) round(x, 1),
                     breaks = seq(-7, 5, by = 1))+
  theme_minimal()# Set breaks at regular intervals


##
labels_for_map=relative_abundances_map %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange(desc(PC1))

relative_abundances_map %>%
  filter(relative_abundance > 0.05, !(object_annotation_category %in% c("multiple organisms", "darksphere"))) %>%
  ggplot(aes(x = as.factor(PC1), y = relative_abundance, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance",
       fill = "Taxa Group") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c(">2 mm","0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short) ->z
z

 
ggsave(
  filename = here("plots/methods_comparison/zooscan_relative_abundances.pdf"), 
  plot = z,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)

taxa_sel=relative_abundances_map
#PLot maps
# Load California map data
worldmap <- map_data("world")
states <- map_data("state")
ca_df <- subset(states, region == "california")




ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=taxa_sel, aes(x=Longitude, y=Latitude,size=relative_abundance, color=size_fraction), alpha=0.7)+ 
  geom_point(data=taxa_sel, aes(x=Longitude, y=Latitude))+
  #Add point at location of max
  scale_size(range = c(2,10))+
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid Copepod Relative Abundance (Zooscan)")+
  theme_classic()+
  facet_wrap(~size_fraction,ncol = 2)

