#Predict proportions using all pools and subpools

##First Size 1

#A1-A3
fido_input_filt_a1_a3_s1=read.csv(here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")%>% 
  select(contains("A1"))

#B3-B5
fido_input_filt_b3_b5_s1=read.csv(here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")%>% 
  select(contains("B3"))


#C5-C7
fido_input_filt_c5_c7_s1=read.csv(here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")%>% 
  select(contains("C5"))

#All
fido_input_filt_all_s1=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family") 

final_s1_input=full_join(rownames_to_column(fido_input_filt_all_s1, var = "RowName"), 
                   rownames_to_column(fido_input_filt_a1_a3_s1, var = "RowName"), 
                   by = "RowName") %>%
  full_join(.,rownames_to_column(fido_input_filt_b3_b5_s1, var = "RowName"),
            by="RowName") %>%
  full_join(.,rownames_to_column(fido_input_filt_c5_c7_s1, var = "RowName"),
            by="RowName") %>%
  # replace(is.na(.), 0) %>%
  column_to_rownames("RowName")
  
                               







#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(final_s1_input))
colnames(final_s1_input) <- gsub("^X", "", colnames(final_s1_input))

##Next, we need to make sure that the orders are the same between meta_18s and final_s1_input
meta_18s <- meta_18s[match(colnames(final_s1_input), meta_18s$Sample_name),] 

meta_18s=meta_18s%>%
  filter(!is.na(Sample_name))

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s1=final_s1_input%>% as.matrix() 



#Fit pibble model 
##MPN: Please remind me how did you choose the 20? Was it using the log marginal likelihood? If so, that code should probably be included here. Happy to chat about this more.
##MPN: This is assuming the default priors for Theta, upsilon, and Xi. Probably reasonable here, but, may want to look in prior predictive checks
##Basically, would run this. These first few rows are just setting the defaults (which fido auto does in the line you have)
upsilon <- nrow(Y_s1)+3 
Omega <- diag(nrow(Y_s1))
G <- cbind(diag(nrow(Y_s1)-1), -1)
Xi <- (upsilon-nrow(Y_s1))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y_s1)-1, nrow(X))
priors <- pibble(NULL, X, Gamma = 20*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
print(priors)
priors <- to_clr(priors)
summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
##Looks ok, centered at zero
##end of added code

##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
fit <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

#Convert to centered log ratio coordinates
fit_s1 <- to_clr(fit)

#Convert to Proportions
fit_prop_1 <- to_proportions(fit_s1)





############







X.tmp.s1 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s1) <- rownames(X)


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s1 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

#Create a dataframe to fill with Cycle 0 proportions thru the loop
final_data_s1 <- data.frame()


#Here begins the loop
for(s in samples_to_loop$sample){
  #Print sample name as a sanity check
  print(s)
  
  #Set selected sample to 1
  X.tmp.s1[s,] <-1
  
  
  #
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
  
  ##MPN: NOt exactly sure what you are trying to show with this plot
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
  
  
  
  
  final_data_s1 <- bind_rows(final_data_s1, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s1[s,] <-1
  
}

#BEEP to notify when finished lol
beepr::beep(12)

#Save final data
current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s1,here(paste0("data/predicted_og/sub_pools/predicted_og_18s_",current_date,"_s1_phy_a1_a3.csv")))



##MPN: Are you aggregating over the samples
final_data_s1 %>% 
  ggplot(., aes(fill=coord, y=n_reads, x=as.factor(cycle_num))) +
  geom_bar(position="stack", stat="identity", width=0.5)+
  scale_fill_discrete(name="ASV")+
  labs(x="PCR Cycle Number",y="Relative Abundance")+
  theme_classic()
## It appears that 'other' is highly over-represented


### Maps for Cycle 0 proportions
#Load complete environmental Metadata file
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

