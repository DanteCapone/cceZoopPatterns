#Preparing inputs for PCR bias-mitigation for COI

library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
library(ggpubr)
library(fido)
library(phyloseq)
here()



#Read in the OTU data
#Run 1 (Non pooled data)
asvcoi_run1=read.csv(here("data/fido/ASV_table_coi_run1.csv")) %>%
  select(-X) 
#Run2
asvcoi_run2=read.csv(here("data/fido/ASV_table_coi_run2.csv")) %>%
  select(-X)



#Taxa Tables 
taxa_coi_meta=read.csv(here("data/past/metazooprunedcoi_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() 

#BlAST
taxa_coi_blast=read.csv(here("data/raw_data/BLAST_taxa_class/leray_taxa.csv")) %>%
  distinct(Hash, .keep_all = TRUE)

taxa_coi=taxa_coi_meta %>% 
  left_join(taxa_coi_blast,., by="Hash") %>%
  select(-contains(".x")) %>%
  rename_all(~gsub("\\.y", "", .)) %>% 
  mutate_all(~replace_na(., "other")) %>% 
  mutate(Hash = if_else(Family == "other", "other", Hash)) %>% 
  distinct() %>% 
  column_to_rownames("Hash")


#Format Long
run1_long=asvcoi_run1 %>%
  pivot_longer(cols = 2:ncol(asvcoi_run1), #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )

run2_long=asvcoi_run2%>%
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
fido_coi_s1=all_runs%>%
  dplyr::select(c(contains("All"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 
fido_coi_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("S2"))) %>% 
  filter(rowSums(.) != 0)
fido_coi_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)

###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_coi_s1_otu=fido_coi_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_coi_s2_otu=fido_coi_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_coi_s3_otu=fido_coi_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa table
taxcoi_s1 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s1))
taxcoi_s1=  tax_table(as.matrix(taxcoi_s1))

#S2
taxcoi_s2 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s2_otu))
taxcoi_s2=  tax_table(as.matrix(taxcoi_s2))
#S3
taxcoi_s3 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s3_otu))
taxcoi_s3=  tax_table(as.matrix(taxcoi_s3))

#Metadata
metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_coi_s1_phy=phyloseq(fido_coi_s1_otu,taxcoi_s1)
fido_coi_s2_phy=phyloseq(fido_coi_s2_otu,taxcoi_s2)
fido_coi_s3_phy=phyloseq(fido_coi_s3_otu,taxcoi_s3)



#PHYLOSEQ
#Agglomerate at the genus level

#S1
fido_coi_s1_phy=phyloseq(fido_coi_s1_otu,taxcoi_s1, metadata)
fido_coi_s1_genus=tax_glom(fido_coi_s1_phy, taxrank = "Genus")

#Check column sums 
colSums(fido_coi_s1_genus_otu)[1:5]
colSums(otu_table(fido_coi_s1_genus))[1:5]

#Make inputs for filtering
fido_coi_s1_genus_otu=otu_table(fido_coi_s1_genus) %>% as.data.frame()
fido_coi_s1_genus_taxa=tax_table(fido_coi_s1_genus) %>% as.data.frame()

#==S2
fido_coi_s2_phy=phyloseq(fido_coi_s2_otu,taxcoi_s2, metadata)
fido_coi_s2_genus=tax_glom(fido_coi_s2_phy, taxrank = "Genus")

#Make inputs for filtering
fido_coi_s2_genus_otu=otu_table(fido_coi_s2_genus) %>% as.data.frame()
fido_coi_s2_genus_taxa=tax_table(fido_coi_s2_genus) %>% as.data.frame()

#==s3
fido_coi_s3_phy=phyloseq(fido_coi_s3_otu,taxcoi_s3, metadata)
fido_coi_s3_genus=tax_glom(fido_coi_s3_phy, taxrank = "Genus")

#Make inputs for filtering
fido_coi_s3_genus_otu=otu_table(fido_coi_s3_genus) %>% as.data.frame()
fido_coi_s3_genus_taxa=tax_table(fido_coi_s3_genus) %>% as.data.frame()


#Save aglomerated family taxa file
taxcoi_genus=rbind(fido_coi_s1_genus_taxa,fido_coi_s2_genus_taxa,fido_coi_s3_genus_taxa) %>%
  unique()

# Create a new row with "other" in all columns
other <- data.frame(lapply(taxcoi_genus, function(x) "other"))

# Bind the new row to the original dataframe
taxcoi_genus <- bind_rows(taxcoi_genus,other)

write.csv(taxcoi_genus,here("data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv"))




#Visualize ECDF
fido_coi_s1_genus_otu[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_coi_s1_genus_otu=fido_coi_s1_genus_otu %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_coi_s1_genus_otu$rowsum, thresh_val)

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_filt <- fido_coi_s1_genus_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_coi_s1_genus_otu %>%
  anti_join(fido_coi_s1_genus_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_coi_s1_final <- bind_rows(fido_taxa_filt, other)


# Remove the identified rows (excluding 'other')
fido_coi_s1_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_coi %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Family,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_coi_s1_save_taxa_phy


#Save
write.csv(fido_coi_s1_save_taxa_phy,here("data/fido/phy/fido_coi_s1_ecdf_taxa_phy.csv"))
write.csv(fido_coi_s1_save_hash_phy,here("data/fido/phy/fido_coi_s1_ecdf_hash_phy.csv"))



### ========== S2 ===============
#s2
fido_coi_s2_genus_otu[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)
###end of added code

fido_coi_s2_genus_otu=fido_coi_s2_genus_otu %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_filt <- fido_coi_s2_genus_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_coi_s2_genus_otu %>%
  anti_join(fido_coi_s2_genus_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_coi_s2_final <- bind_rows(fido_taxa_filt, other)


# Remove the identified rows (excluding 'other')
fido_coi_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_coi %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  select(-rowsum,-Phylum,-Class,-Family,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_coi_s2_save_taxa_phy

#Hash
# select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_coi_s2_save_hash_phy

#Save
write.csv(fido_coi_s2_save_taxa_phy,here("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv"))
write.csv(fido_coi_s2_save_hash_phy,here("data/fido/phy/fido_coi_s2_ecdf_hash_phy.csv"))





### =============== S3 =========
#s3
fido_coi_s3_genus_otu[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)
###end of added code

fido_coi_s3_genus_otu=fido_coi_s3_genus_otu %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_filt <- fido_coi_s3_genus_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_coi_s3_genus_otu %>%
  anti_join(fido_coi_s3_genus_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_coi_s3_final <- bind_rows(fido_taxa_filt, other)




# Remove the identified rows (excluding 'other')
fido_coi_s3_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_coi %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  select(-rowsum,-Phylum,-Class,-Family,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_coi_s3_save_taxa_phy

#Hash
# select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_coi_s3_save_hash_phy

#Save
write.csv(fido_coi_s3_save_taxa_phy,here("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv"))
write.csv(fido_coi_s3_save_hash_phy,here("data/fido/phy/fido_coi_s3_ecdf_hash_phy.csv"))

