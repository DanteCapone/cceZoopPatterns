# From these runs, we didn't get data from the 20C or various pools for COI runs, but it worked for 18s
# try first with 18s

library (tidyverse)
library (here)
library(ggpubr)
library(fido)
library(phyloseq)
here()



###18S
#Read in the OTU data
#Run 1 (Non pooled data)
asv18s_run1=read.csv(here("data/past/","ASV_table_18s_run1.csv")) %>%
  select(-X) 
#Run2
asv18s_run2=read.csv(here("data/past/","ASV_table_18s_run2.csv")) %>%
  select(-X)



#Taxa Tables 
taxa_18s_meta=read.csv(here("data/past/metazoopruned18s_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() 

#BlAST
taxa_18s_blast=read.csv(here("data/raw_data/BLAST_taxa_class/zhang_taxa.csv")) %>%
distinct(Hash, .keep_all = TRUE)

taxa_18s=taxa_18s_meta %>% 
  left_join(taxa_18s_blast,., by="Hash") %>%
  select(-contains(".x")) %>%
  rename_all(~gsub("\\.y", "", .)) %>% 
  mutate_all(~replace_na(., "other")) %>% 
  mutate(Hash = if_else(Family == "other", "other", Hash)) %>% 
  distinct() %>% 
  column_to_rownames("Hash")



filter# 2) Merging and manipulation (updated 8/24/2023 to create a new 18S input for fido where
# I don't average technical replicates)

#Format Long
run1_long=asv18s_run1 %>%
  pivot_longer(cols = 2:ncol(asv18s_run1), #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )

run2_long=asv18s_run2%>%
  pivot_longer(cols = 2:ncol(.),  #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )



all_runs=bind_rows(run1_long,run2_long) %>%
  pivot_wider(names_from = Sample_ID, values_from = Nreads)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  column_to_rownames("Hash")

#Replace X
colnames(all_runs) <- gsub("^X", "", colnames(all_runs))


#Separate out by size
#S1
fido_18s_s1=all_runs%>%
  dplyr::select(c(contains("All"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 
fido_18s_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("S2"))) %>% 
  filter(rowSums(.) != 0)
fido_18s_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)

# MPN: Just to clarify, what is the difference between the ".1" and ".2" samples?
# E.g., C1.T7.H9.1 versus C1.T7.H9.2?


###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_18s_s1_otu=fido_18s_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu=fido_18s_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu=fido_18s_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa tables
tax18s_s1 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>% 
  mutate(Family = ifelse(Family=="", 'other', Family))
tax18s_s1=  tax_table(as.matrix(tax18s_s1))

#S2
tax18s_s2 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>% 
  mutate(Family = ifelse(Family=="", 'other', Family))
tax18s_s2=  tax_table(as.matrix(tax18s_s2))
#S3
tax18s_s3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>% 
  mutate(Family = ifelse(Family=="", 'other', Family))
tax18s_s3=  tax_table(as.matrix(tax18s_s3))






#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1)
fido_18s_s2_phy=phyloseq(fido_18s_s2_otu,tax18s_s2)
fido_18s_s3_phy=phyloseq(fido_18s_s3_otu,tax18s_s3)



#PHYLOSEQ
#Agglomerate at the family level

#S1
fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1, metadata)
fido_18s_s1_family=tax_glom(fido_18s_s1_phy, taxrank = "Family")

#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s1))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s1_family))

# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s1_family_otu=otu_table(fido_18s_s1_family) %>% as.data.frame() 
  
fido_18s_s1_family_taxa=tax_table(fido_18s_s1_family) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s1_family_otu <- bind_rows(fido_18s_s1_family_otu, difference)


colSums(fido_18s_s1_otu)[1:5]
colSums(fido_18s_s1_family_otu)[1:5]

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s1_family_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s1_family_taxa) -> fido_18s_s1_family_taxa



#==S2
fido_18s_s2_phy=phyloseq(fido_18s_s2_otu,tax18s_s2, metadata)
fido_18s_s2_family=tax_glom(fido_18s_s2_phy, taxrank = "Family")


#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s2))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s2_family))

# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s2_family_otu=otu_table(fido_18s_s2_family) %>% as.data.frame()
fido_18s_s2_family_taxa=tax_table(fido_18s_s2_family) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s2_family_otu <- bind_rows(fido_18s_s2_family_otu, difference)


colSums(fido_18s_s2_otu)[1:5]
colSums(fido_18s_s2_family_otu)[1:5]


#Make inputs for filtering
fido_18s_s2_family_otu=otu_table(fido_18s_s2_family) %>% as.data.frame()
fido_18s_s2_family_taxa=tax_table(fido_18s_s2_family) %>% as.data.frame()

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s2_family_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s2_family_taxa) -> fido_18s_s2_family_taxa



#==s3
fido_18s_s3_phy=phyloseq(fido_18s_s3_otu,tax18s_s3, metadata)
fido_18s_s3_family=tax_glom(fido_18s_s3_phy, taxrank = "Family")


#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s3))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s3_family))

# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s3_family_otu=otu_table(fido_18s_s3_family) %>% as.data.frame() 
  
fido_18s_s3_family_taxa=tax_table(fido_18s_s3_family) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s3_family_otu <- bind_rows(fido_18s_s3_family_otu, difference)


colSums(fido_18s_s3_otu)[1:5]
colSums(fido_18s_s3_family_otu)[1:5]



#Make inputs for filtering
fido_18s_s3_family_otu=otu_table(fido_18s_s3_family) %>% as.data.frame()
fido_18s_s3_family_taxa=tax_table(fido_18s_s3_family) %>% as.data.frame()

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s3_family_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s3_family_taxa) -> fido_18s_s3_family_taxa

#Save aglomerated family taxa file
tax18s_family=rbind(fido_18s_s1_family_taxa,fido_18s_s2_family_taxa,fido_18s_s3_family_taxa) %>%
  unique() %>%
  mutate(
    Kingdom = if_else(Family == 'other', 'other', Kingdom),
    Phylum = if_else(Family == 'other', 'other', Phylum),
    Subphylum = if_else(Family == 'other', 'other', Subphylum),
    Class = if_else(Family == 'other', 'other', Class),
    Subclass = if_else(Family == 'other', 'other', Subclass),
    Superorder = if_else(Family == 'other', 'other', Superorder),
    Order = if_else(Family == 'other', 'other', Order),
  ) %>% 
  select(-Species, -Genus)
write.csv(tax18s_family,here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))


## MPN: I think there is a bug here. Why does colSums(fido_18s_s1_final) and colSums(fido_18s_s1) not match?

## ==== S1 ====
# Separate rows based appearance in the calibration samples
<<<<<<< Updated upstream
# MPN: The filter seems strict (not that it matters looking at the data). You are only allowing 0 counts in two samples or less (for the calibration samples)?
fido_taxa_filt <- fido_18s_s1_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
=======
fido_taxa_filt <- fido_18s_s1_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")
>>>>>>> Stashed changes
other <- fido_18s_s1_family_otu %>%
  anti_join(fido_18s_s1_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(Hash = "other")

# Combine data
fido_18s_s1_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")


colSums(fido_18s_s1_final)[1:5]


#Join with taxa file
fido_18s_s1_save_family_phy <- fido_18s_s1_final %>%
  rownames_to_column("Hash") %>%
  left_join(fido_18s_s1_family_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


#Save
write.csv(fido_18s_s1_save_family_phy,here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"))



## ==== s2 ====
fido_taxa_filt <- fido_18s_s2_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")
other <- fido_18s_s2_family_otu %>%
  anti_join(fido_18s_s2_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(Hash = "other")

# Combine data
fido_18s_s2_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

#Join with taxa file
fido_18s_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s2_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s2_save_family_phy

#Save
write.csv(fido_18s_s2_save_family_phy,here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv"))





## ==== s3 ====
# Separate rows based appearance in the calibration samples
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
# Separate rows based appearance in the calibration samples
fido_taxa_filt <- fido_18s_s3_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")
other <- fido_18s_s3_family_otu %>%
  anti_join(fido_18s_s3_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(Hash = "other")

# Combine data
fido_18s_s3_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

#Join with taxa file
fido_18s_s3_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s3_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s3_save_family_phy

#Save
write.csv(fido_18s_s3_save_family_phy,here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv"))

