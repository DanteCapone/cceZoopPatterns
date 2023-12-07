#Size based fido model
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(ggpubr)
library(fido)
library(stringr)
library(here)
library(gridExtra)


fido_input_filt=read.csv(file.path("data/fido/fido_18s_s1_ecdf_spp_hash.csv"), header=TRUE, check.names = FALSE, row.names = 1)
  
  #Metadata
  meta_18s=read.csv(file.path("data/fido","meta_18s_unaveraged_s1.csv"), header=TRUE) %>%
    select(-c(X)) %>%
    filter(Sample_name %in% colnames(fido_input_filt))
  colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))
  
  ##MPN: Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
  meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]
  
  #Model matrix
  ##MPN: To be clear, this will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
  X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
  
  Y_s1=fido_input_filt%>% as.matrix() 
  
  ## MPN: Cleaning the names slightly to make it easier to read
  i <- 1:nrow(Y_s1)
  rownames(Y_s1) <- sub("NA", "", rownames(Y_s1))
  rownames(Y_s1) <- paste0("seq_", i, "_", rownames(Y_s1))
  rownames(Y_s1) <- sub("^(.*\\..*\\..{5}).*", "\\1",(rownames(Y_s1)))
  
  fit <- pibble(Y_s1, X, gamma = 20*diag(nrow(X)), n_samples = 10000)
  
  # ,Convert to centered log ratio coordinates
  fit_s1 <- to_clr(fit)

############
###Proportions

#Predict at cycle 0
fit_prop_1 <- to_proportions(fit_s1)
# predicted_s1 <- predict(fit_prop_1, newdata=X.tmp.s1, summary=TRUE) %>% 



#Select sample to predict on
#Sample select
sample_sel="sample_numC1.T7.H9_S1"



X.tmp.s1 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s1) <- rownames(X)


#Samples to loop thru
X.tmp.s1 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

final_data_s1 <- data.frame()

for(s in samples_to_loop$sample){
  print(s)
  X.tmp.s1[s,] <-1
  
  predicted_s1 <- predict(fit_prop_1, newdata=X.tmp.s1, summary=TRUE) %>% 
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("0.2-0.5mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s1 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>%
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", ""))))%>%
    rownames_to_column("coord") %>%
    pivot_longer(cols = c(-coord),
                 names_to = "replicate",
                 values_to = "n_reads") %>%
    mutate(size=rep("0.2-0.5mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s1,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
  

  
  ###
  # taxa_list=rownames(Y_s1)
  # taxa_sel=taxa_list[1:10]
  # focus.coord <- paste0("clr_", taxa_sel) 
  # focus.covariate <- rownames(X.tmp.s1)[which(grepl("sample_num", rownames(X.tmp.s1)))]
  # ##  
  # predicted_s1 %>% filter(coord %in% focus.coord) %>% 
  #   ggplot(aes(x=cycle_num)) +
  #   geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey") +
  #   geom_line(aes(y=mean)) +
  #   geom_point(data=tidy_calibration %>% filter(coord %in% focus.coord), aes(y=val)) +
  #   facet_grid(coord~.) +
  #   theme_bw() +
  #   theme(strip.text.y=element_text(angle=0)) +
  #   ylab("CLR Coordinates") ##Again, would be to proportions, Not CLR coordinates :)
  # 
  
  final_data_s1 <- bind_rows(final_data_s1, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s1[s,] <-1
  
}


beepr::beep(12)

write.csv(final_data_s1,here("data/predicted_og/predicted_og_18s_11_3_2023_s1.csv"))


##Barplots of predicted C0 proportions
# taxa_list=predicted_s1 %>%
#   arrange(desc(n_reads)) %>%
#   select(coord)
# taxa_sel=taxa_list[1:10,]
# focus.coord <- taxa_sel
# 
# predicted_s1 %>% filter(coord %in% focus.coord) %>%
#   ggplot(., aes(fill=coord, y=mean, x=as.factor(cycle_num))) + 
#   geom_bar(position="stack", stat="identity", width=0.5)+
#   scale_fill_discrete(name="ASV")+
#   labs(x="PCR Cycle Number",y="Relative Abundance")+
#   facet_wrap(~size, nrow=3)+
#   theme_classic()



### Maps for OG proportions
#Metadata
metazoo_meta=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>% 
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot")

metazoo_meta_map=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>% 
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  mutate(offshore_onshore=metazoo_meta$offshore_onshore)%>%
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
  dplyr::select(Latitude,Longitude) %>%
  rownames_to_column("sample_id")


#Now make a dataframe for mapping and add lat/long
map_pcr_18s_s1=final_data_s1 %>% filter(cycle_num==0)%>%
  mutate(sample_id = str_extract(replicate, "(C|CT)\\d+\\.T\\d+\\.H\\d+_S\\d+")) %>%
  left_join(metazoo_meta_map, by="sample_id")%>%
  filter(grepl("Calanoida", coord, ignore.case = TRUE))



#Load metadata
# Load California map data
worldmap <- map_data("world")
states <- map_data("state")
ca_df <- subset(states, region == "california")


p1=ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=map_pcr_18s_s1, aes(x=Longitude, y=Latitude,size=n_reads, color=size), alpha=0.7)+ 
  # geom_point(data=taxa_sel_all, aes(x=object_lon, y=object_lat,size=concentraion, color=size_fraction, alpha=0/5))+ 
  #Add point at location of max
  scale_size(range = c(2,12))+
  geom_point(data=map_pcr_18s_s1, aes(x=Longitude, y=Latitude))+ 
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid ASV Relative Read Abundance 
       (18S PCR bias-mitigated ASV's)")+
  theme_classic()
p1








###############
############### Let's repeat for other sizes now

############First 0.5-1############
fido_input_filt=read.csv(file.path("data/fido/fido_18s_s2_ecdf_spp_hash.csv"), header=TRUE, check.names = FALSE, row.names = 1)

#Metadata
meta_18s=read.csv(file.path("data/fido","meta_18s_unaveraged_s2.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##MPN: Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
##MPN: To be clear, this will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))

Y_s2=fido_input_filt%>% as.matrix() 

## MPN: Cleaning the names slightly to make it easier to read
i <- 1:nrow(Y_s2)
rownames(Y_s2) <- sub("NA", "", rownames(Y_s2))
rownames(Y_s2) <- paste0("seq_", i, "_", rownames(Y_s2))
rownames(Y_s2) <- sub("^(.*\\..*\\..{5}).*", "\\1",(rownames(Y_s2)))

fit <- pibble(Y_s2, X, gamma = 20*diag(nrow(X)), n_samples = 10000)

# ,Convert to centered log ratio coordinates
fit_s2 <- to_clr(fit)

############
###Proportions

#Predict at cycle 0
fit_prop_2 <- to_proportions(fit_s2)
# predicted_s2 <- predict(fit_prop_1, newdata=X.tmp.s2, summary=TRUE) %>% 


X.tmp.s2 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s2) <- rownames(X)


#Samples to loop thru
X.tmp.s2 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

final_data_s2 <- data.frame()

for(s in samples_to_loop$sample){
  print(s)
  X.tmp.s2[s,] <-1
  
  predicted_s2 <- predict(fit_prop_2, newdata=X.tmp.s2, summary=TRUE) %>% 
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("0.5-1mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s2 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>%
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", ""))))%>%
    rownames_to_column("coord") %>%
    pivot_longer(cols = c(-coord),
                 names_to = "replicate",
                 values_to = "n_reads") %>%
    mutate(size=rep("0.5-1mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s2,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
 
  final_data_s2 <- bind_rows(final_data_s2, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s2[s,] <-1
  
}


beepr::beep(4)

write.csv(final_data_s2,here("data/predicted_og/predicted_og_18s_11_3_2023_s2.csv"))

### Maps for OG proportions

#Now make a dataframe for mapping and add lat/long
map_pcr_18s_s2=final_data_s2 %>% filter(cycle_num==0)%>%
  mutate(sample_id = str_extract(replicate, "(C|CT)\\d+\\.T\\d+\\.H\\d+_S\\d+")) %>%
  left_join(metazoo_meta_map, by="sample_id")%>%
  filter(grepl("Calanoida", coord, ignore.case = TRUE))



#Load metadata
# Load California map data
worldmap <- map_data("world")
states <- map_data("state")
ca_df <- subset(states, region == "california")


p2=ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=map_pcr_18s_s2, aes(x=Longitude, y=Latitude,size=n_reads, color=size), alpha=0.7)+ 
  # geom_point(data=taxa_sel_all, aes(x=object_lon, y=object_lat,size=concentraion, color=size_fraction, alpha=0/5))+ 
  #Add point at location of max
  scale_size(range = c(2,12))+
  geom_point(data=map_pcr_18s_s2, aes(x=Longitude, y=Latitude))+ 
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid ASV Relative Read Abundance 
       (18S PCR bias-mitigated ASV's)")+
  theme_classic()
p2




#########
##### 1-2mm####
fido_input_filt=read.csv(file.path("data/fido/fido_18s_s3_ecdf_spp_hash.csv"), header=TRUE, check.names = FALSE, row.names = 1)

#Metadata
meta_18s=read.csv(file.path("data/fido","meta_18s_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##MPN: Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
##MPN: To be clear, this will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))

Y_s3=fido_input_filt%>% as.matrix() 

## MPN: Cleaning the names slightly to make it easier to read
i <- 1:nrow(Y_s3)
rownames(Y_s3) <- sub("NA", "", rownames(Y_s3))
rownames(Y_s3) <- paste0("seq_", i, "_", rownames(Y_s3))
rownames(Y_s3) <- sub("^(.*\\..*\\..{5}).*", "\\1",(rownames(Y_s3)))

fit <- pibble(Y_s3, X, gamma = 20*diag(nrow(X)), n_samples = 10000)

# ,Convert to centered log ratio coordinates
fit_s3 <- to_clr(fit)
# plot(fit_s3, par="Lambda", focus.cov="cycle_num")
X.tmp.s3 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s3) <- rownames(X)


###Proportions

#Predict at cycle 0
fit_prop_3 <- to_proportions(fit_s3)
# predicted_s3 <- predict(fit_prop_1, newdata=X.tmp.s3, summary=TRUE) %>% 




#Make X.tmp for loop
X.tmp.s3 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s3) <- rownames(X)


#Samples to loop thru
X.tmp.s3 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

final_data_s3 <- data.frame()

for(s in samples_to_loop$sample){
  print(s)
  X.tmp.s3[s,] <-1
  
  predicted_s3 <- predict(fit_prop_3, newdata=X.tmp.s3, summary=TRUE) %>% 
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("1-2mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s3 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>%
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", ""))))%>%
    rownames_to_column("coord") %>%
    pivot_longer(cols = c(-coord),
                 names_to = "replicate",
                 values_to = "n_reads") %>%
    mutate(size=rep("1-2mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s3,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
  
  
  
  final_data_s3 <- bind_rows(final_data_s3, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s3[s,] <-1
  
}


beepr::beep(12)

write.csv(final_data_s3,here("data/predicted_og/predicted_og_18s_11_3_2023_s3.csv"))


### Maps for OG proportions
#Now make a dataframe for mapping and add lat/long

final_data_s3=read.csv(here("data/predicted_og/predicted_og_18s_11_3_2023_s3.csv")) %>%
  select(-X.1,X)


map_pcr_18s_s3=final_data_s3 %>% filter(cycle_num==0)%>%
  mutate(sample_id = str_extract(replicate, "(C|CT)\\d+\\.T\\d+\\.H\\d+_S\\d+")) %>%
  left_join(metazoo_meta_map, by="sample_id") %>%
  filter(grepl("Calanoida", coord, ignore.case = TRUE))



#Load metadata
# Load California map data
worldmap <- map_data("world")
states <- map_data("state")
ca_df <- subset(states, region == "california")


p3=ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=map_pcr_18s_s3, aes(x=Longitude, y=Latitude,size=n_reads),color="green", alpha=0.7)+ 
  # geom_point(data=taxa_sel_all, aes(x=object_lon, y=object_lat,size=concentraion, color=size_fraction, alpha=0/5))+ 
  #Add point at location of max
  scale_size(range = c(2,12))+
  geom_point(data=map_pcr_18s_s3, aes(x=Longitude, y=Latitude))+ 
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid ASV Relative Read Abundance 
       (18S PCR bias-mitigated ASV's)")+
  theme_classic()
p3


#Final plot

#Join all 3
final_data_all_sizes=rbind(final_data_s1,final_data_s2,final_data_s3)

map_pcr_18s=final_data_all_sizes %>% filter(cycle_num==0)%>%
  mutate(sample_id = str_extract(replicate, "(C|CT)\\d+\\.T\\d+\\.H\\d+_S\\d+")) %>%
  left_join(metazoo_meta_map, by="sample_id") %>%
  filter(grepl("Calanoida", coord, ignore.case = TRUE)) %>%
  group_by(Latitude, Longitude, sample_id) %>%
  summarize(
    n_reads = sum(n_reads),
    size = toString(unique(size))
  )
  
ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(map_id=region), col = "white", fill = "gray50") +
  geom_point(data=map_pcr_18s, aes(x=Longitude, y=Latitude,size=n_reads,color=size), alpha=0.7)+ 
  # geom_point(data=taxa_sel_all, aes(x=object_lon, y=object_lat,size=concentraion, color=size_fraction, alpha=0/5))+ 
  #Add point at location of max
  scale_size(range = c(2,12))+
  geom_point(data=map_pcr_18s, aes(x=Longitude, y=Latitude))+ 
  coord_fixed(xlim = c(-134, -119.0),  ylim = c(34, 38), ratio = 1.3)+
  scale_x_continuous(breaks = seq(-118,-132, by = -2))+
  xlab("Latitude")+
  ylab("Longitude")+
  labs(title="Calanoid ASV PCR Bias-Corrected Relative Read Abundance (18S)")+
  facet_wrap(~size, ncol=2)+
  theme_classic()

