#Zooscan analysis final
library(here)
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




#Compute relative abudances 
relative_abundances=zooscan_processed %>%
  group_by(sample_id,size_fraction,object_annotation_category) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count),
         relative_abundance = count / total)

##Lat and lon are funky so add back in from metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_")))


# Standardize key columns (some samples are incorerctly notated to compare with Zooscan)
metadata <- metadata %>%
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_")))

# Merge dataframes: add metadata to relative abudance data
relative_abundances_map <- relative_abundances %>%
  left_join(metadata, by="sample_id")




####PLOTTING

###Zooscan relative abundances vs. pc1

ggplot(relative_abundances_map, aes(x=))



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

