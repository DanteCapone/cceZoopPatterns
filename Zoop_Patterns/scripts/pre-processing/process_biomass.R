#Add depth
depths=read.csv(here("PCR_bias_correction/data/physical_environmental_data/sample_depths.csv")) %>%
  select(-X) %>%
  mutate(Sample_ID_short=Sample_ID)


#Volume filtered (add to metadata)
volume_filtered=read.csv(here("PCR_bias_correction/data/raw_data/biomass/p2107_bt_volume_filtered.csv"))

#Dryweights
dryweights=read.csv("PCR_bias_correction/data/raw_data/biomass/dryweights_forzoopmetab.csv") %>%
  mutate(biomass_dry=8/3*biomass_dry) %>% 
  left_join(.,volume_filtered, by = c("Sample_ID_short"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>% 
  #Replace missing tow depths with 210
  mutate(depth = replace(depth, which(is.na(depth)), 210)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID) 

write.csv(dryweights,here("Zoop_Patterns/data/zoop_other/biomass_processed.csv"))
