## Data Descriptions:
The data here include 3 types: amplicon sequence variant (ASV) tables, taxonomy tables, and pre-processed dataframes for fido input. There were sample plate PCR runs: the first run was environmental samples run in duplicate while the second contained a mix of sample replicates and pooled samples for the variable PCR experiment. Sample IDs are as follows:
EX) C1-T7-H9_S1.1 = Cycle#-Tow #- Haul #_Size class.Replicate
Cycle, tow and haul refer to plankton net sampling site and number accoridng to cruise operations. Size class is 1 of 3 size classes (0.2-0.5, 0.5-1, or 1-2mm) we separated plankton into to study ecolological/life history characteristics. 

For Pooled Samples we had 4 different pools that are identified by the samples that were aggregated into that pool. The "A1-A3", "B3-B5", "C5-C7" refer the the 96-well plate grid coordinates from the PCR plate that contains the sequence of samples aggregated. The "All" pool contains an aliquot from all samples.

Some of the files 

### Formatted/Processed Files
These files have been formatted from the raw output ASV tables such that:
-Row names have been creaetd by combining the Hash code with the best available taxonomic resolution
-Hashes that don't have a count >1 in more than 30% of samples have been removed 
-The sample replicates have been averaged 
-Hashes which don't occur in all of the pooled samples are combined into an "Other" category (these have been subsetted to the correct pool so I separated the samples into the "A1","B3" and a hash must occur in every sample and replicate of the variaible PCR (20,24,28)...)

*fido_18s_all.csv*: File for the pool with all samples
*fido_18s_a1.csv*: File for the A1-A3 pool
*fido_18s_b3.csv*: File for the B3-B5 pool
*fido_18s_c5.csv*: File for the C5-C7 pool

The smaller pools contain the following samples

A1-A3: *"C1-T7-H9_S1", "C1-T7-H9_S2", "C1-T7-H9_S3", "C1-T8-H10_S1", "C1-T8-H10_S2", "C1-T8-H10_S3", "C2-T8-H18_S1"
          , "C2-T8-H18_S2", "C2-T8-H18_S3", "C2-T9-H19_S1", "C2-T9-H19_S2", "C2-T9-H19_S3", "C3-T6-H25_S1", 
          "C3-T6-H25_S2", "C3-T6-H25_S3", "C3-T7-H26_S1", "C3-T7-H26_S2"*

B3-B5: *"C3-T7-H26_S3", "CT1-T1-H28_S1", "CT1-T1-H28_S2", "CT1-T1-H28_S3"
          , "CT1-T2-H29_S1", "CT1-T2-H29_S2", "CT1-T2-H29_S3", "CT1-T3-H30_S1", "CT1-T3-H30_S2", "CT1-T3-H30_S3", "CT1-T4-H31_S1", "CT1-T4-H31_S2"
          , "CT1-T4-H31_S3", "CT1-T5-H32_S1", "CT1-T5-H32_S2", "CT1-T5-H32_S3", "CT1-T6-H33_S1"*

C5-C7: *"CT1-T6-H33_S2", "CT1-T6-H33_S3", "CT1-T7-H34_S1", "CT1-T7-H34_S2", "CT1-T7-H34_S3",
          "CT1-T8-H35_S1", "CT1-T8-H35_S2", "CT1-T8-H35_S3", "CT2-T1-H36_S1", "CT2-T1-H36_S2", "CT2-T1-H36_S3",
          "CT2-T4-H39_S1", "CT2-T4-H39_S2", "CT2-T4-H39_S3", "CT2-T8-H43_S1", "CT2-T8-H43_S2", "CT2-T8-H43_S3"*

### Fido metadata
Files with the metadata as required by the fido tutorial

*meta_18s_c5.csv*: Contains sample names, PCR cycle number, etc

### Raw ASV files
These are files straight from the pipeline output. The data is a matrix where row is associated with a particular ASV code/Hash and the row is associated with a sample replicate. There are 2 runs for each gene marker
*ASV_table_18s_run1.csv*: This file contains the raw datamatrix for the 1st run from 18S gene marker. This is just environmental data and not the variaible PCR experiment pools

*ASV_table_18s_run2.csv*: This file contains the raw datamatrix for the 2nd run from 18S gene marker. This is contains the variable PCR pool experiment samples as well as some environmental replicates.




### Taxonomy Files
*metazoopruned18s_tax.csv*: 18S rRNA marker metabarcoding taxonomy file classified using metazoogene database
*metazooprunedcoi_tax.csv*: COI marker metabarcoding taxonomy file classified using metazoogene database


# Main Data (Environmental Samples, not Variable PCR Experiment)

## Phyloseq Bio Data

### 18S
These are datafiles associaetd with the 18S primer, which is the primer that worked properly for the variable PCR experiment. There are 2 files that can be used with the _phyloseq_ package to make a phyloseq object when combined with the _env_metadata_impute_phyloseq_6.9.2023.csv_ file in the physical environmental data folder
1. OTU table
2. Taxa table


### COI
These are datafiles associaetd with the COI primer, which provides better species level identification when compared to 18S. There are 2 files that can be used with the _phyloseq_ package to make a phyloseq object when combined with the _env_metadata_impute_phyloseq_6.9.2023.csv_ file in the physical environmental data folder
1. OTU table
2. Taxa table



## Physical Environmental Data 
There are several different versions of the physical environmental metadata file in this folder
*Note* there are only 17 unique sampling sites, but each site has 3 samples (though some failed in which case there are 2. If using unique sites only, I index the rows using "Sample_ID_short" but if using the complete data to compare to biological analysis/phyloseq I use "Sample_ID".


1. _env_metadata_all_samples_8.16.2023.csv_*: This file has metadata assocaited with each unique biological sample in the dataset, however given there are 3 samples per physical site (one for each zooplankton size) there is redundacy. This contains all the continuous environmental variaibles used in the PCA 
2. _env_metadata_impute_phyloseq_6.9.2023.csv*_: This file is used for the phyloseq objects and contains additional non-continuous variables and results from PCA/clustering including: Sizefractionmm, max_size, cycle, Sample_ID_short, Sample_ID_dot, PC1, clust_group, and offshore_onshore
3. _env_metadata_just_sites_8.16.2023.csv**_: This is the same as 1, but only contains unique values for each sampling site
4. _env_pca_dist_for_clustering_8.8.2023.csv**_: This contains the actual PCA values from the PCA analysis and can be fed directly into the hierarchical clustering

*Contains data for biological samples (51)
**Contains only unique data for sampling sites (17)
