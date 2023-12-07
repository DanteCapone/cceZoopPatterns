#Pie plot on map of top taxa for 18s/coi onshore/offshore


#Import packages
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(maps)
library(dplyr)
library(tidyverse)
library(here)
library(phyloseq)

#Functions for plotting
ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)

source(("scripts/treemap_funs_Capone.R"))
source("scripts/phyloseq_mapping_funs.R")


#Read in the data
#COI
leray_metazoo_otucoi=read.csv(here("data/Phyloseq Objects/COI pseq/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash")%>%
  dplyr::select(where(~ !is.na(.[[1]])))
column_to_rownames("Sample_ID_dot")
leray_metazoo_taxa=read.csv(here("data/phyloseq_objects_eDNAindex/coi_taxa_table_eDNA_metazoogene.csv")) %>% column_to_rownames("X")


#18s
zhan_otu=read.csv(here("data/Phyloseq Objects/18s pseq/metazoopruned18s_otu.csv"), header=TRUE)%>%
  select(where(~ !any(is.na(.)))) %>%
  column_to_rownames("Hash")
zhan_taxa=read.csv(here("data/Phyloseq Objects/18s pseq/metazoopruned18s_tax.csv ")) %>%
  column_to_rownames("Hash")

#Metadata
metazoo_meta=read.csv(here("data/CURRENT_WORKING_Metadata/env_metadata_impute_phyloseq_6.9.2023.csv"))%>% 
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot")
metazoo_meta_map=read.csv(here("data/CURRENT_WORKING_Metadata/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>% 
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  mutate(offshore_onshore=metazoo_meta$offshore_onshore)




#Make phyloseq objects
OTU = otu_table(as.matrix(leray_metazoo_otucoi), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(leray_metazoo_taxa))
meta=sample_data(metazoo_meta_map)
Phy_norm_coi <- phyloseq_normalize_median(phyloseq(OTU, TAX, meta)) %>%
  #Transform to proportions
  transform_sample_counts(., function(x) x / sum(x) )%>%
  phyloseq_transform_to_long(.)


#18s
OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(metazoo_meta_map)
Phy_norm_18s <- phyloseq_normalize_median(phyloseq(OTU, TAX, meta)) %>%
  #Transform to proportions
  transform_sample_counts(., function(x) x / sum(x) ) %>%
  phyloseq_transform_to_long(.)







# Load California map data
worldmap <- map_data("world")
states <- map_data("state")
ca_df <- subset(states, region == "california")




# Calculate category counts and sort the levels

#Select phyloseq object
phy_sel=Phy_norm_coi

#Group by calanoida and sample

#Load taxa from PCR bias mitigated data for comparability

##COI
coi_taxa_pcr=read.csv(here("data/taxa_lists/coi_pcrmitigated_taxa.csv")) %>% select(-"X")
phy_sel %>% filter(asv_code %in% coi_taxa_pcr$Hash) %>%
  filter(Order=="Calanoida") %>%
  group_by(asv_code,file_code,Latitude,Longitude,Sizefractionmm,cycle) %>%
  group_by(file_code,Latitude,Longitude,Sizefractionmm,cycle) %>%
  summarize(n_reads=sum(n_reads))%>%
  filter_all(all_vars(!is.na(.)))%>%
  mutate(Sizefractionmm= ifelse(Sizefractionmm == "2-Jan", "1-2",Sizefractionmm))->calanoid_pcr_coi

##18S
zhan_taxa_pcr=read.csv(here("data/taxa_lists/zhan_pcrmitigated_taxa.csv")) %>% select(-"X")
phy_sel=Phy_norm_18s
phy_sel %>% filter(asv_code %in% zhan_taxa_pcr$Hash) %>%
  filter(Order=="Calanoida") %>%
  group_by(asv_code,file_code,Latitude,Longitude,Sizefractionmm,cycle) %>%
  group_by(file_code,Latitude,Longitude,Sizefractionmm,cycle) %>%
  summarize(n_reads=sum(n_reads))%>%
  filter_all(all_vars(!is.na(.)))%>%
  mutate(Sizefractionmm= ifelse(Sizefractionmm == "2-Jan", "1-2",Sizefractionmm))->calanoid_pcr_18s

#Save dataframes
# write.csv(calanoid_pcr_18s,here("data/methods_compare/calanoid_pcr_18.csv"),row.names = FALSE)
# write.csv(calanoid_pcr_coi,here("data/methods_compare/calanoid_pcr_coi.csv"),row.names = FALSE)


##PLOT

##COI
ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=calanoid_pcr, aes(x=Longitude, y=Latitude,size=n_reads, color=Sizefractionmm), alpha=0.7)+ 
  # geom_point(data=taxa_sel_all, aes(x=object_lon, y=object_lat,size=concentraion, color=size_fraction, alpha=0/5))+ 
  #Add point at location of max
  scale_size(range = c(2,12))+
  geom_point(data=calanoid_pcr, aes(x=Longitude, y=Latitude))+ 
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid ASV Relative Read Abundance of 
       (Calanus pacificus & Ctenocalanus Vanus ASV's)")+
  theme_classic()+
  facet_wrap(~Sizefractionmm,ncol = 2)


##18s
##COI
ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=calanoid_pcr_18s, aes(x=Longitude, y=Latitude,size=n_reads, color=Sizefractionmm), alpha=0.7)+ 
  scale_size_continuous(range = c(2,12), breaks = c(0.2,0.4,0.6), labels = c(0.2, 0.4, 0.6)) +
  # geom_point(data=taxa_sel_all, aes(x=object_lon, y=object_lat,size=concentraion, color=size_fraction, alpha=0/5))+ 
  #Add point at location of max
  # scale_size(range = c(2,12))+
  geom_point(data=calanoid_pcr_18s, aes(x=Longitude, y=Latitude))+ 
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid Raw ASV Relative Read Abundance (18S)")+
  theme_classic()+
  facet_wrap(~Sizefractionmm,ncol = 2)












phy_sel$asv_code <- factor(phy_sel$asv_code, levels = sorted_levels)


# Calculate the counts of the specific column
counts <- phy_sel %>%
  dplyr::count(asv_code)




#Filter to asv's with more than n counts in all 48 samples
phy_sel=phy_sel %>% group_by(asv_code) %>%
  filter(n() > 24) %>%
  filter(Species != "NA") %>%
  ungroup() 


#Add variabe corresponding to histogram of occurence
phy_sel_counts=phy_sel%>%
  group_by(asv_code) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(rank = rownames(.))

#Add ranks to the original df
phy_sel= phy_sel%>%
  left_join(phy_sel_counts, by = "asv_code") 




#Histogram
ggplot(phy_sel_counts, aes(x = asv_code, y = count, fill = asv_code)) +
  geom_bar(stat = "identity") +
  labs(title = "Histogram of Counts", x = "Variable", y = "Count") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate the percentages for each category in each region

asv_1=phy_sel %>% filter(rank==1)%>%
  #group by cycle, size fraction and hash
  group_by(cycle,asv_code,max_size) %>% 
  #Sum each hash occurence
  mutate(sum_asv=sum(n_reads)) %>%
  ungroup() %>%
  #sum each hash at each cycle
  group_by(Latitude,max_size) %>%
  mutate(nreads_cycle=sum(n_reads)) %>%
  mutate(percent = sum_asv/nreads_cycle)%>%
  mutate(eDNA.Index=nreads_cycle)





#Map ASV Normalized 
map_asv_norm(phy_sel,1)

map_asv_norm(phy_sel,2)
map_asv_norm(phy_sel,3)
map_asv_norm(phy_sel,4)
map_asv_norm(phy_sel,5)
map_asv_norm(phy_sel,6)


#Raw
map_asv_raw(phy_raw_map,1)
map_asv_raw(phy_raw_map,2)
map_asv_raw(phy_raw_map,3)
ggsave()
map_asv_raw(phy_raw_map,4)
map_asv_raw(phy_raw_map,5)
map_asv_raw(phy_raw_map,6)



#### 7/9/2023
#Look at 5 selected asv's

top3_coi=read.csv(here("data/maps/coi_top3.csv")) %>%
  mutate(asv_code=Hash)
top3_18s=read.csv(here("data/maps/18s_top30.csv"))%>%
  mutate(asv_code=Hash)

#Filter to top 3
phy_sel_top3_coi=Phy_norm_coi %>% 
  filter(asv_code %in% top3_coi$Hash)%>%
  left_join(top3_coi %>% select(rank_off_on,rank,offshore_onshore, asv_code), by="asv_code") %>%
  rename(offshore_onshore.y="offshore_onshore_taxa")

#eDNA index
phy_sel_top3_coi_eDNA=Phy_edna_coi %>% 
  filter(asv_code %in% top3_coi$Hash)%>%
  left_join(top3_coi %>% select(rank_off_on,rank,offshore_onshore, asv_code), by="asv_code") %>%
  rename(offshore_onshore.y="offshore_onshore_taxa")


phy_sel_top3_18S=Phy_norm_18s %>% 
  filter(asv_code %in% top3_18s$Hash)%>%
  left_join(top3_18s %>% select(rank_off_on,rank,offshore_onshore, asv_code), by="asv_code") %>%
  rename(offshore_onshore.y="offshore_onshore_taxa") %>%
  #Add in family if spp is missing for 18s
  mutate(Species = ifelse(Species == "", NA, Species)) %>%
  mutate(Species = ifelse(is.na(Species), paste0("Unidentified ",Family), Species))


asv_codes_coi=unique(phy_sel_top3_coi$asv_code)
asv_codes_18s=unique(phy_sel_top3_18s$asv_code)



#COI
#Map
p1=map_asv_off_on(phy_sel_top3_coi,"1_offshore","COI")
ggsave("figures/Maps/top3_onshore_offshore_taxa/coi_offshore_rank1_map_7102023.png", plot = p1, width = 11, height = 8.4, units = "in")
