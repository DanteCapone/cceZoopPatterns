# From these runs, we didn't get data from the 20C or various pools for COI runs, but it worked for 18s
# try first with 18s

library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
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
taxa_18s=read.csv(here("data/past/metazoopruned18s_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>%
  column_to_rownames("Hash")

#BlAST
taxa_18s=read.csv(here("data/raw_data/BLAST_taxa_class/zhang_taxa.csv")) %>% 
  distinct(Hash, .keep_all = TRUE)%>% column_to_rownames("Hash")


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



###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_18s_s1_otu=fido_18s_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu=fido_18s_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu=fido_18s_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa tables
tax18s_s1 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1))
tax18s_s1=  tax_table(as.matrix(tax18s_s1))

#S2
tax18s_s2 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu))
tax18s_s2=  tax_table(as.matrix(tax18s_s2))
#S3
tax18s_s3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu))
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

#Make inputs for filtering
fido_18s_s1_family_otu=otu_table(fido_18s_s1_family) %>% as.data.frame()
fido_18s_s1_family_taxa=tax_table(fido_18s_s1_family) %>% as.data.frame() 

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
  unique()
write.csv(tax18s_family,here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))




## ==== S1 ====
# Separate rows based appearance in the calibration samples
fido_taxa_filt <- fido_18s_s1_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s1_family_otu %>%
  anti_join(fido_18s_s1_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")

# Combine data
fido_18s_s1_final <- bind_rows(fido_taxa_filt, other)

#Join with taxa file
fido_18s_s1_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s1_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  #Hash
  select(-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash)->fido_18s_s1_save_family_phy

#Save
write.csv(fido_18s_s1_save_family_phy,here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"))



## ==== s2 ====
# Separate rows based appearance in the calibration samples
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_filt <- fido_18s_s2_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s2_family_otu %>%
  anti_join(fido_18s_s2_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")

# Combine data
fido_18s_s2_final <- bind_rows(fido_taxa_filt, other)

#Join with taxa file
fido_18s_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s2_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  #Hash
  select(-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash)->fido_18s_s2_save_family_phy

#Save
write.csv(fido_18s_s2_save_family_phy,here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv"))





## ==== s3 ====
# Separate rows based appearance in the calibration samples
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_filt <- fido_18s_s3_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s3_family_otu %>%
  anti_join(fido_18s_s3_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")

# Combine data
fido_18s_s3_final <- bind_rows(fido_taxa_filt, other)

#Join with taxa file
fido_18s_s3_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s3_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  #Hash
  select(-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash)->fido_18s_s3_save_family_phy

#Save
write.csv(fido_18s_s3_save_family_phy,here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv"))

