#Zooscan functions

readEcotaxa <- function(data) {
  #data: A .tsv file exported from Ecotaxa
  #select relevant columns to calculate data
  #selected columns for database upload
  data$sample <- data$sample_id
  data$Haul <- data$sample_id
  data$Region <- "California Current"

  data$Detail_Location <- data$sample_id
  data$Comment <- ""
  data$process_particle_pixel_size_mm <- 0.0106
  
  
  #Caclulate volume filtered to get concentration
  data =data %>%
    mutate(sample_conc=acq_sub_part/sample_tot_vol) %>%
    mutate(cycle= str_extract(object_id, "^[^_-]+"))
  
  #Select relevant columns
  data_select <- data %>% dplyr::select(., sample_ship, sample_program, sample_id, Haul, Region, Detail_Location, Comment, 
                                                  object_date, object_time, object_lat, object_lon, sample_bottomdepth, object_depth_min,
                                                  object_depth_max, object_annotation_category, object_annotation_hierarchy, object_annotation_person_name,
                                                  sample, object_id, sample_id, sample_tot_vol, acq_sub_part,object_feret, 
                                                  object_area, object_major, object_minor, object_area_exc, process_particle_pixel_size_mm, acq_max_mesh,
                                                  sample_conc,cycle, object_annotation_status,acq_id)
  #Filter: only living and validated
  data_select<-data_select %>% 
    filter(object_annotation_status=="validated")%>%
    filter(!str_detect(object_annotation_hierarchy, regex("not-living", ignore_case = TRUE)))
  
  
  #Convert to mm for calculating ESD
  data_select$area_mm2  <- data_select$object_area * (data_select$process_particle_pixel_size_mm**2) 
  
  data_select$major_mm  <- data_select$object_major * data_select$process_particle_pixel_size_mm
  
  data_select$minor_mm  <- data_select$object_minor * data_select$process_particle_pixel_size_mm
  
  data_select$area_exc_mm2  <- data_select$object_area_exc * (data_select$process_particle_pixel_size_mm**2) 
  
  data_select$area_majmin_mm2  <- pi * data_select$major_mm/2 * data_select$minor_mm/2
  
  data_select$esd_mm  <- 2 * (sqrt(data_select$area_mm2/pi))
  
  data_select$esd_exc_mm  <- 2 * (sqrt(data_select$area_exc_mm2/pi))
  
  data_select$esd_maj_min_mm <- 2 * (sqrt(data_select$area_majmin_mm2/pi))
  
  print(unique(data_select$object_annotation_category))
  
  #Extract size fraction
  data_final=data_select%>%
    mutate(size_fraction=as.factor(acq_max_mesh))
  
  #Add size group as a category
  data_final %>%
    mutate(size_fraction = case_when(
      esd_mm >= 0.2 & esd_mm < 0.5  ~ '0.2-0.5',
      esd_mm >= 0.5 & esd_mm < 1    ~ '0.5-1',
      esd_mm >= 1   & esd_mm < 2    ~ '1-2',
      esd_mm > 2                          ~ '>2',
      TRUE                                      ~ 'Other'
    ))-> data_final
  
  
 return(data_final)

  
}
# 
# 
# #Load in the data
# 
# 
# biovolume <- read.table(here("data/Zooscan/ecotaxa_export_8030_20230823_2344.tsv"),header=TRUE, sep="\t", encoding="latin1")
# 
# #select relevant columns to calculate biovolume
# #selected columns for database upload
# biovolume$sample <- biovolume$sample_id
# biovolume$Haul <- biovolume$sample_id
# biovolume$Region <- "California Current"
# biovolume$Detail_Location <- biovolume$sample_id
# biovolume$Comment <- ""
# biovolume$process_particle_pixel_size_mm <- 0.0106
# 
# 
# #Caclulate volume filtered to get concentration
# biovolume =biovolume %>%
#   mutate(sample_conc=acq_sub_part/sample_tot_vol) %>%
#   mutate(cycle= str_extract(object_id, "^[^_-]+"))
# 
# 
# biovolume_select <- biovolume %>% dplyr::select(., sample_ship, sample_program, sample_id, Haul, Region, Detail_Location, Comment, 
#                                                 object_date, object_time, object_lat, object_lon, sample_bottomdepth, object_depth_min,
#                                                 object_depth_max, object_annotation_category, object_annotation_hierarchy, object_annotation_person_name,
#                                                 sample, object_id, sample_id, sample_tot_vol, acq_sub_part,object_feret, 
#                                                 object_area, object_major, object_minor, object_area_exc, process_particle_pixel_size_mm, acq_max_mesh,
#                                                 sample_conc,cycle, object_annotation_status,acq_id)
# 
# 
# 
# 
# #delete all non-plankton categories (adjust according to dataset)
# # biovolume_select<-biovolume_select[!(biovolume_select$object_annotation_category=="bubble" | 
# #                                        biovolume_select$object_annotation_category=="fiber<detritus"| 
# #                                        biovolume_select$object_annotation_category=="Ellobiopsidae"|
# #                                        biovolume_select$object_annotation_category=="multiple<plastic"|
# #                                        biovolume_select$object_annotation_category=="multiple<other"|
# #                                        biovolume_select$object_annotation_category=="detritus" | 
# #                                        biovolume_select$object_annotation_category=="egg<Acartia sinjiensis" | 
# #                                        biovolume_select$object_annotation_category=="artefact"| 
# #                                        biovolume_select$object_annotation_category=="antenna<Crustacea" | 
# #                                        biovolume_select$object_annotation_category=="leg<Crustacea"| 
# #                                        biovolume_select$object_annotation_category=="dead<Copepoda"| 
# #                                        biovolume_select$object_annotation_category=="Ostracoda X"| 
# #                                        biovolume_select$object_annotation_category=="egg sac<egg"| 
# #                                        biovolume_select$object_annotation_category=="feces" |
# #                                        biovolume_select$object_annotation_category=="part<Copepoda" |
# #                                        biovolume_select$object_annotation_category=="Foraminifera"),]
# 
# 
# 
# #Extract size fraction
# biovolume_final=biovolume_select%>%
#   mutate(size_fraction=as.factor(acq_max_mesh))
# 
# 
# #Convert to mm
# biovolume_final$area_mm2  <- biovolume_final$object_area * (biovolume_final$process_particle_pixel_size_mm**2) 
# 
# biovolume_final$major_mm  <- biovolume_final$object_major * biovolume_final$process_particle_pixel_size_mm
# 
# biovolume_final$minor_mm  <- biovolume_final$object_minor * biovolume_final$process_particle_pixel_size_mm
# 
# biovolume_final$area_exc_mm2  <- biovolume_final$object_area_exc * (biovolume_final$process_particle_pixel_size_mm**2) 
# 
# biovolume_final$area_majmin_mm2  <- pi * biovolume_final$major_mm/2 * biovolume_final$minor_mm/2
# 
# biovolume_final$esd_mm  <- 2 * (sqrt(biovolume_final$area_mm2/pi))
# 
# biovolume_final$esd_exc_mm  <- 2 * (sqrt(biovolume_final$area_exc_mm2/pi))
# 
# biovolume_final$esd_maj_min_mm <- 2 * (sqrt(biovolume_final$area_majmin_mm2/pi))
# 
# print(unique(biovolume_final$object_annotation_category))
# 
# 
# 
# 
# biovolume_final %>%
#   #Add size group
#   mutate(size_fraction = case_when(
#     esd_mm >= 0.2 & esd_mm < 0.5  ~ '0.2-0.5',
#     esd_mm >= 0.5 & esd_mm < 1    ~ '0.5-1',
#     esd_mm >= 1   & esd_mm < 2    ~ '1-2',
#     esd_mm > 2                          ~ '>2',
#     TRUE                                      ~ 'Other'
#   ))-> biovolume_final
# 
# #Validated
# biovolume_validated=biovolume_final %>%
#   filter(object_annotation_status=="validated")%>%
#   filter(!str_detect(object_annotation_hierarchy, regex("not-living", ignore_case = TRUE)))
# 
# relative_abundances=biovolume_validated %>%
#   group_by(sample_id,size_fraction,object_annotation_category) %>%
#   summarise(count = n()) %>%
#   mutate(total = sum(count),
#          relative_abundance = count / total)
# 
# ## Clean up some issues
# relative_abundances=relative_abundances %>%
#   mutate(sample_id = str_replace_all(sample_id, "-", "_"),
#          sample_id = ifelse(sample_id == "c2_t1_h36", "ct2_t1_h36", sample_id),
#          sample_id = ifelse(sample_id == "ct1_t8_h10", "c1_t8_h10", sample_id),
#          sample_id = ifelse(sample_id == "ct2_t9_h19", "c2_t9_h19", sample_id),
#          sample_id = ifelse(sample_id == "c3_bt6_h25", "c3_t6_h25", sample_id))
# 
# ##Lat and lon are funky so add back in from metadata
# metadata=read.csv(here("data/CURRENT_WORKING_Metadata/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))
# 
# # Standardize key columns
# metadata <- metadata %>%
#   mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
#   dplyr::select(sample_id,Latitude,Longitude)
# 
# # Merge dataframes
# relative_abundances_map <- relative_abundances %>%
#   left_join(metadata, by="sample_id")
# 
# 
# #Select a taxa
# taxa_pick="Calanoida"
# taxa_sel=relative_abundances_map %>%
#   filter(object_annotation_category==taxa_pick) %>%
#   filter(size_fraction != ">2")
# 
# # taxa_sel_valid=taxa_sel %>%
# #   filter(object_annotation_status=="validated") %>%
# #   group_by(object_lat,object_lon,sample_id, size_fraction,sample_conc) %>%
# #   summarize(count = n()) %>%
# #   mutate(concentraion=sample_conc*count) 
# #   
# # taxa_sel_all=taxa_sel %>%
# #   group_by(object_lat,object_lon,sample_id, size_fraction,sample_conc) %>%
# #   summarize(count = n()) %>%
# #   mutate(concentraion=sample_conc*count) 
# 
# 
