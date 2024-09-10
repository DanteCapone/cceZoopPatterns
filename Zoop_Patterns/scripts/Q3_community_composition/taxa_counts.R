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



# #Taxa Tables 
# taxa_18s=read.csv(here("data/past/metazoopruned18s_tax.csv"))%>%
#   mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
#   group_by(Hash) %>%
#   filter(rank(desc(non_na_count)) == 1) %>%
#   select(-non_na_count) %>%
#   ungroup() %>%
#   column_to_rownames("Hash")

#BlAST
taxa_18s=read.csv(here("data/raw_data/BLAST_taxa_class/zhang_taxa.csv")) %>% 
  distinct(Hash, .keep_all = TRUE)%>% column_to_rownames("Hash") 

taxa_18s=tax_table(as.matrix(taxa_18s))

taxa_18s

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
  column_to_rownames("Hash")%>%
  dplyr::select(c(contains("All"),contains("S"))) 

#Replace X
colnames(all_runs) <- gsub("^X", "", colnames(all_runs))

all_runs = all_runs %>% 
  otu_table(taxa_are_rows = TRUE)



#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)


blast_18s_phy=phy=phyloseq(all_runs,taxa_18s)

# Assuming physeq is your phyloseq object

# Extract taxonomy table
taxa_table <- tax_table(blast_18s_phy) %>% as.data.frame()

# Initialize an empty list to store the results
result_list <- list()

# Loop through each taxonomic level
for (level in colnames(taxa_table)) {
  # Count unique values at the current taxonomic level
  taxa_count <- table(taxa_table[, level])
  
  # Append the result to the list
  result_list[[level]] <- length(unique(names(taxa_count)))
}

# Combine the results into a data frame
result_df <- data.frame(Taxa_Level = names(result_list), Unique_Categories = unlist(result_list))

# Display the final data frame
print(result_df)



# COI ---------------------------------------------------------------------

# ,Read in the data
# ,Run 1 (Non pooled data)
asvcoi_run1=read.csv(here("data","fido","ASV_table_coi_run1.csv")) %>% 
  dplyr::select(-X)


# ,Run2
asvcoi_run2=read.csv(here("data","fido","ASV_table_coi_run2.csv"))%>% 
  dplyr::select(-X)


# ,Taxa Tables 
taxa_coi=read.csv(here("data/raw_data/BLAST_taxa_class/leray_taxa.csv")) %>%
  column_to_rownames("Hash")

taxa_coi=tax_table(as.matrix(taxa_coi))
  



# 2) Merging and manipulation (updated 8/24/2023 to create a new coi input for fido where
# I don't average technical replicates)
# First need to average technical replicates
# To do this i need to format long
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

all_runs = all_runs %>% 
  otu_table(taxa_are_rows = TRUE)



#Metadata
metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)


blast_coi_phy=phy=phyloseq(all_runs,taxa_coi)

# Assuming physeq is your phyloseq object

# Extract taxonomy table
taxa_table <- tax_table(blast_coi_phy) %>% as.data.frame()

# Initialize an empty list to store the results
result_list <- list()

# Loop through each taxonomic level
for (level in colnames(taxa_table)) {
  # Count unique values at the current taxonomic level
  taxa_count <- table(taxa_table[, level])
  
  # Append the result to the list
  result_list[[level]] <- length(unique(names(taxa_count)))
}

# Combine the results into a data frame
result_df <- data.frame(Taxa_Level = names(result_list), Unique_Categories = unlist(result_list))

# Display the final data frame
print(result_df)


