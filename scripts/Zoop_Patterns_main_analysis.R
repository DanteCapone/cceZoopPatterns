# Zooplankton CCE Barcoding Paper 1 Analysis Script

########### Load the libraries
library(tidyverse)
library(factoextra)
library(gridExtra)
library(ggfortify)
library(ggrepel)
library(here)
library(missMDA)
library(factoextra)
library(NbClust)
library(corrplot)
library(matrixStats)
library(ggpubr)
library(vegan)
library(corrplot)
library(ggnewscale)
library(phyloseq)
library(rerddap)
library(RColorBrewer)

#Project path
project_path="."

########## Load the helper scripts
source(here("scripts/helpful_functions/treemap_funs_Capone.R"))
source(here("scripts/helpful_functions/phyloseq_mapping_funs.R")
source(here("scripts/helpful_functions/general_helper_functions.R")




# #Load the data ----------------------------------------------------------

#Metadata
#1) load in metadata and select desired rows 
env_metadata_raw<-read.csv(here("data/pre_processing/metadata03-28-2023.csv")) 

#Part one: Read in the physical environmental data

env_metadata<-read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv")) %>% dplyr::select(-c("X"))%>%
  column_to_rownames("Sample_ID_dot") %>%
  mutate(PC1=PC1*-1) %>% 
  #Add day night for plotting
  mutate(day_night = ifelse(day_night_0_1 == 0, "day", "night"))
  


#Physloseq

#Load in the phyloseq data and format to a table

#COI raw reads
coi_otu=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
coi_meta=env_metadata 

coi_taxa=read.csv(here("data/taxa_files/blast_metazoo_coi.csv")) %>% column_to_rownames("Hash")


#Merge by site
coi_meta_all=coi_meta %>% group_by(Sample_ID_short) %>%
  summarize_all(median) %>% 
  column_to_rownames("Sample_ID_short")

OTU = otu_table(as.matrix(coi_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(coi_meta)
meta$cycle=meta$cycle %>% as.factor()
Phy_coi_raw <- phyloseq(OTU, TAX, meta)


#Merge by sample site
Phy_merged_coi_raw <- merge_samples(Phy_coi_raw,c("Sample_ID_short"))


#Issues with metadata...create a new phyloseq object with merged metadata
Phy_merged_coi_raw=phyloseq(otu_table(Phy_merged_coi_raw),tax_table(Phy_merged_coi_raw),sample_data(coi_metazoo_meta_all))

#18s
zhan_otu=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
zhan_taxa=read.csv(here("data/taxa_files/blast_metazoo_18s.csv")) %>% column_to_rownames("Hash")
zhan_meta=env_metadata


#Merge by site
zhan_meta_all=zhan_meta %>% group_by(Sample_ID_short) %>%
  summarize_all(median) %>% 
  column_to_rownames("Sample_ID_short")

OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(zhan_meta)
meta$cycle=meta$cycle %>% as.factor()
Phy_zhan_raw <- phyloseq(OTU, TAX, meta)

#Merge by sample site
Phy_merged_zhan_raw <- merge_samples(Phy_zhan_raw,c("Sample_ID_short"))


#Issues with metadata...create a new phyloseq object with merged metadata
Phy_merged_zhan_raw=phyloseq(otu_table(Phy_merged_zhan_raw),tax_table(Phy_merged_zhan_raw),sample_data(zhan_meta_all))








# Main Analysis ###########################################


######### Q1 Physical Analysis: #Script for Map, EDA visualization, PCA, and clustering -------------------


# #1) Code for SS Chlorophyll Map -----------------------------------------




# Set color palette (similar to the color range in your image)
color_breaks <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)  # Adjusted breaks

# Set color palette
my_palette <- brewer.pal(11, "RdYlBu")  # 11 colors for a diverging color scale

#Reverse the color palette
reversed_palette <- rev(my_palette)



# Coordinates for plotting
coordinates <- data.frame(
  Latitude = c(
    33.606036, 36.512291, 36.511306, 36.169252, 36.163166,
    36.142257, 36.140615, 36.181659, 36.186195, 36.186618,
    36.113802, 36.098667, 36.116328, 36.266852, 36.293449,
    36.356773, 36.380195, 36.407122, 36.431205, 36.429056,
    36.460569, 34.751460, 34.713717, 34.670477, 34.631567,
    34.589691, 34.560463, 34.528097, 34.503944, 34.960436,
    35.049822, 35.023270, 35.095029, 35.116252, 35.157629,
    35.203037, 35.178398, 35.563177, 35.429799, 35.282737,
    35.131118, 34.990777, 34.826846, 34.695320, 34.542449
  ),
  Longitude = c(
    -119.323358, -122.278297, -122.276371, -122.039465, -122.061647,
    -122.130000, -122.223000, -122.337500, -122.414167, -122.466667,
    -122.495500, -122.446000, -122.367500, -122.643667, -122.709500,
    -122.827500, -122.899333, -122.993500, -123.065667, -123.091167,
    -123.225667, -130.548333, -130.539833, -130.525833, -130.518333,
    -130.510315, -130.477410, -130.459688, -130.423656, -128.575522,
    -127.572833, -126.635167, -125.523500, -124.721998, -123.885993,
    -123.105333, -122.300500, -121.606667, -121.504167, -121.385833,
    -121.260667, -121.161833, -121.025167, -120.926833, -120.828500
  )
)

# Set site coordinates and time for chlorophyll data extraction
chlStartDate <- "2021-07-13"
chlEndDate <- "2021-08-13"

# Adjusted longitudes to 0-360 range
chlCoordsLon <- c(-132, -116)  # Adjusted longitude to match your plot
chlCoordsLat <- c(31, 40)      # Adjusted latitude to match your plot

# Set dataset source for chlorophyll data
chlSource <- info("erdMBchla1day")

# Get chlorophyll data
chlorophyll <- griddap(chlSource,
                       time = c(chlStartDate, chlEndDate),
                       longitude = chlCoordsLon+360,
                       latitude = chlCoordsLat,
                       fields = "chlorophyll",
                       fmt = "csv")

# Log-transform the chlorophyll data
chlorophyll$log_chlorophyll <- log(chlorophyll$chlorophyll)

# Monthly composite for California region
chlorophyll$month <- format(as.Date(chlorophyll$time), "%Y-%m")

monthly_composite_chl <- chlorophyll %>%
  group_by(latitude, longitude, month) %>%
  summarize(mean_log_chlorophyll = mean(log_chlorophyll, na.rm = TRUE)) %>%
  ungroup()

# Load map data for California
california_map <- map_data("state", region = "california")

# Create the map using ggplot
chl_map=ggplot() +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "black") +
  geom_tile(data = monthly_composite_chl, aes(x = longitude-360, y = latitude, fill = mean_log_chlorophyll), na.rm = TRUE) +
  scale_fill_gradientn(colours = reversed_palette, values = scales::rescale(log(color_breaks)), 
                       limits = log(c(0.05, 15)), breaks = log(color_breaks), labels = color_breaks,
                       name = "Log Chlorophyll") +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Monthly Composite log Sea Surface Chlorophyll") +
  coord_fixed(1.3, xlim = c(-131, -118), ylim = c(32, 38)) +
  theme(
    axis.text = element_text(size = 14),        # Increase axis tick size
    axis.title = element_text(size = 16, face = "bold"),  # Increase axis label size
    plot.title = element_text(size = 18, face = "bold"),  # Increase title size
    legend.text = element_text(size = 12),      # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

#Add stations



# Adding station points, add day night as the shape

#Load colors
my_palette=custom_pallete_all()

chl_map_p2107=chl_map+
  new_scale_fill() +
  geom_point(data = env_metadata, aes(x = Longitude, y = Latitude, shape = day_night, fill = Sample_ID_short), size = 10) +
  scale_shape_manual(values = c("day" = 21, "night"=23), name= "Day/Night") +
  scale_fill_manual(values = my_palette, name="Station ID")+
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 2)) +  # Set fill legend to 2 columns
  theme(legend.title = element_text(face = "bold"))  # Make sure legend titles are styled consistently


chl_map_p2107
  

# Save
#PNG & PDF Save
ggsave(
  filename = here(project_path,"figures/Q1_physical_analysis/sampling_map_chla.png"),
  plot = chl_map_p2107,
  width = 12,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here(project_path,"figures/Q1_physical_analysis/sampling_map_chla.pdf"),
  plot = chl_map_p2107,
  width = 12,  # Width in inches
  height = 8  # Height in inches
)




# Counting taxa in each level -----------------------------------------




# 18S ---------------------------------------------------------------------


#Read in the OTU data
#Run 1 (Non pooled data)
asv18s_run1=read.csv(here(project_path,"data/raw_reads/","ASV_table_18s_run1.csv")) %>%
  select(-X) 
#Run2
asv18s_run2=read.csv(here(project_path,"data/raw_reads/","ASV_table_18s_run2.csv")) %>%
  select(-X)

#BlAST
taxa_18s=read.csv(here(project_path,"data/raw_data/BLAST_taxa_class/zhang_taxa.csv")) %>% 
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
meta18s=read.csv(here(project_path,"data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)


blast_18s_phy=phy=phyloseq(all_runs,taxa_18s)


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
asvcoi_run1=read.csv(here(project_path,"data/raw_reads/","ASV_table_coi_run1.csv")) 


# Run2
asvcoi_run2=read.csv(here(project_path,"data/raw_reads/","ASV_table_coi_run2.csv"))

all_runs=bind_rows(asvcoi_run1,asvcoi_run2) %>%
  pivot_wider(names_from = Sample_name, values_from = nReads)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  column_to_rownames("Hash")

all_runs = all_runs %>% 
  otu_table(taxa_are_rows = TRUE)


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







#Metadata
metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)


blast_coi_phy=phy=phyloseq(all_runs,taxa_coi)


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




#2) EDA and PCA


# Making the metadata -----------------------------------------------------

#Various metadata frames for each analysis
env_metadata_long= env_metadata_raw %>% dplyr::select(c(Sample_ID,max_size,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                    PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                    PAR_1_depth)) 
env_metadata_phy= env_metadata_raw %>% dplyr::select(c(Sample_ID_dot,Sample_ID,max_size,Sizefractionmm,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                    PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                    PAR_1_depth))
env_metadata_phy_map= env_metadata_raw %>% dplyr::select(c(Sample_ID_dot,Sample_ID,max_size,Sizefractionmm,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                   PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                   PAR_1_depth,LONGITUDE,LATITUDE))

#If I want to just look at sampling site not all sizes since the data are repeated for each size
#Make nighttime PAR equal to day time
env_metadata_sel =env_metadata_long %>%
  group_by(cycle) %>%
  mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>%
  filter(max_size==0.5) %>%
  dplyr::select(-max_size)%>% 
  bind_rows() %>% ungroup() %>%
  dplyr::select(-c("Sample_ID","cycle","PAR_1_depth"))%>%column_to_rownames("Sample_ID_short") 

#####BIOLOGICAL METADATA Make an imputed dataset for phyloseq analysis

# env_metadata_phy_sel =env_metadata_phy %>%
#   group_by(cycle) %>%
#   mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>% 
#   bind_rows() %>% ungroup() %>%
#   dplyr::select(-c("Sample_ID_short","cycle","PAR_1_depth","max_size","Sizefractionmm","Sample_ID_dot"))%>%column_to_rownames("Sample_ID")
# 
# #For cluster analysis use short DF 
env_metadata_phy_sel =env_metadata_phy %>%
  group_by(cycle) %>%
  mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>%
  bind_rows() %>% ungroup() %>%
  dplyr::select(-c("Sample_ID","cycle","PAR_1_depth","max_size","Sizefractionmm","Sample_ID_dot"))%>%
  distinct()%>%column_to_rownames("Sample_ID_short")




# IMPUTATION --------------------------------------------------------------

## Imputation, impute the missing nutrient data
estim_ncpPCA(env_metadata_sel)
metadata_impute <- imputePCA(env_metadata_sel,ncp=1)

## Imputation on phyloseq metadata frames
estim_ncpPCA(env_metadata_phy_sel)
metadata_impute <- imputePCA(env_metadata_phy_sel,ncp=1)

env_metadata_phy_add_back =env_metadata_raw %>%
  filter(max_size==0.5)
metadata_impute_df_phy=metadata_impute$completeObs %>% data.frame() %>% 
  #Add back in categoraical vars
  mutate(Sizefractionmm=env_metadata_phy_add_back$Sizefractionmm) %>%
  mutate(max_size=env_metadata_phy_add_back$max_size) %>%
  mutate(cycle=env_metadata_phy_add_back$cycle) %>%
  mutate(Sample_ID_short=env_metadata_phy_add_back$Sample_ID_short) %>%
  mutate(Sample_ID_dot=env_metadata_phy_add_back$Sample_ID_dot)   












# Correlation plot --------------------------------------------------------

here()
env_metadata_corr=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X"))  %>% 
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size))%>% 
  dplyr::select(-Sample_ID_short)

##Env Correlation matrix
corr_res=env_metadata_corr %>% 
  cor(.)

# Compute p-values using correlation matrix
p_values <- cor.mtest(env_metadata_corr)$p %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(cols = -variable, names_to = "variable2", values_to = "p.value")%>%
  # Adjust p-values using Benjamini-Hochberg correction
  mutate(p.adj= p.adjust(p.value, method = "BH") )

p_vals_adj=p_values %>%
  group_by(variable, variable2) %>%
  summarise(p.adj = mean(p.adj, na.rm = TRUE)) %>%
  ungroup() %>% # Ensure to ungroup the data after summarizing 
  pivot_wider(names_from = variable, values_from = p.adj) %>%
  column_to_rownames("variable2") %>%
  as.matrix()




# Corr using corrplot --------------------------------------------------------------------
dev.off()
# Sort the row names and column names to ensure alignment
corr_res <- corr_res[order(rownames(corr_res)), order(colnames(corr_res))]
p_vals_adj <- p_vals_adj[order(rownames(p_vals_adj)), order(colnames(p_vals_adj))]


#Make correlation plot: blue is positive red is negative
corr_plot=corrplot(corr_res,p.mat=p_vals_adj, type = 'lower', order = 'FPC', tl.col = 'black',
                   cl.ratio = 0.2, tl.srt = 45)
corr_plot

#PNG & PDF Save
ggsave(
  filename = here("plots/Q1_physical_analysis/corr_plot_p_adj.png"),
  plot = corr_plot,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/Q1_physical_analysis/corr_plot_p_adj.pdf"),
  plot = corr_plot,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)




# PCA ---------------------------------------------------------------------

pca_in=metadata_impute$completeObs %>% as.data.frame() 
env_pca=prcomp(pca_in, scale=TRUE)
# env_pca=PCA(pca_in)
summary(env_pca)
env_pca_df=env_pca$x %>% as.data.frame()



#New plot
env_pca_df=env_pca$x %>% as.data.frame()



#New plot

##Contribution of vars to pcs
#Which PCs contribute significantly
explained_var <- env_pca$sdev^2 / sum(env_pca$sdev^2)
df <- data.frame(component = 1:length(explained_var), explained_var = explained_var)
fviz_eig(env_pca, addlabels = TRUE)

var<-get_pca_var(env_pca)
pc1_p<-fviz_contrib(env_pca,"var",axes = 1) # default angle=45?
plot(pc1_p,main = "Variables percentage contribution of PC1")


var<-get_pca_var(env_pca)
pc2_p<-fviz_contrib(env_pca,"var",axes = 2) # default angle=45?
plot(pc2_p,main = "Variables percentage contribution of PC2")
grid.arrange(pc1_p,pc2_p)


#Using ggfortify

#Biplot
fviz_pca_biplot(env_pca, 
                pointsize = 3,  # size of data points
                repel = TRUE)
fviz_pca_var(env_pca, col.var = "cos2",
             gradient.cols = c("blue" ,"purple","orange","red"),
             repel = TRUE)+
  theme(
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())


## Ellipses

plot1=autoplot(env_pca, data=metadata_impute_df_phy,label=TRUE,colour="cycle", label.size=4, loadings=TRUE,size=3, loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size =5, scale = 1, repel=TRUE)+theme_classic()
plot1



#Clustering using pvclust. Cluster the PCA dataframe using hierarchical clustering
# 
library(pvclust)
set.seed(123)
env_pca_df=env_pca$x %>% as.data.frame()


#Determine the Optimal cluster #
# Silhouette method
fviz_nbclust(pca_in, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method") 
#Result: 2 clusters is optimal

#hierarchical clustering 
res.pv <- pvclust(t(env_pca_df), method.dist="cor",method.hclust="average", nboot = 10000)

#Plot and save figures
plot(res.pv, hang = -1, cex = 1.2,xlab="Sample ID", main="Environmental PCA Clusterings \n with p-values (%)")
pvrect(res.pv, alpha = 0.95, lwd=2)

#PNG & PDF Save
ggsave(
  filename = here("plots/Q1_physical_analysis/pca_clusters.png"),
  plot = corr_plot,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/Q1_physical_analysis/pca_clusters.pdf"),
  plot = corr_plot,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)





res.pv.phys=res.pv %>% as.data.frame(.)




#Add Clusters to metadata
clusters <- pvpick(res.pv)
clusters
onshore=array(unlist(clusters$clusters[2]),dim=c(8,1))
offshore=array(unlist(clusters$clusters[1]),dim=c(9,1))
onshore



#Add clustering variable to dataframe
metadata_impute_df_phy$clust_group=rep(1,length(metadata_impute_df_phy$Sample_ID_short))
metadata_impute_df_phy$clust_group[metadata_impute_df_phy$Sample_ID_short %in% offshore]=2
metadata_impute_df_phy$offshore_onshore[metadata_impute_df_phy$clust_group==2]="offshore"
metadata_impute_df_phy$offshore_onshore[metadata_impute_df_phy$clust_group==1]="onshore"



#Add PC1 to the metadata file
metadata_impute_df_phy=left_join(metadata_impute_df_phy, env_pca_df %>% 
                                   select(Sample_ID_short,PC1), by="Sample_ID_short")






# Q2: Diversity Patterns Shannon vs. PC1 --------------------------------------------------------------------

#Load in the phyloseq data and format to a table

coi_metazoo_meta=env_metadata 

#Merge by site
coi_metazoo_meta_all=coi_metazoo_meta %>% group_by(Sample_ID_short) %>%
  summarize_all(median) %>% 
  column_to_rownames("Sample_ID_short")

OTU = otu_table(as.matrix(coi_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(coi_metazoo_meta)
meta$cycle=meta$cycle %>% as.factor()
Phy_coi_raw <- phyloseq(OTU, TAX, meta)


#Merge by sample site
Phy_merged_coi_raw <- merge_samples(Phy_coi_raw,c("Sample_ID_short"))


#Issues with metadata...create a new phyloseq object with merged metadata
Phy_merged_coi_raw=phyloseq(otu_table(Phy_merged_coi_raw),tax_table(Phy_merged_coi_raw),sample_data(coi_metazoo_meta_all))

#18s
zhan_meta=env_metadata

#Merge by site
zhan_meta_all=zhan_meta %>% group_by(Sample_ID_short) %>%
  summarize_all(median) %>% 
  column_to_rownames("Sample_ID_short")

OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(zhan_meta)
meta$cycle=meta$cycle %>% as.factor()
Phy_zhan_raw <- phyloseq(OTU, TAX, meta)

#Merge by sample site
Phy_merged_zhan_raw <- merge_samples(Phy_zhan_raw,c("Sample_ID_short"))


#Issues with metadata...create a new phyloseq object with merged metadata
Phy_merged_zhan_raw=phyloseq(otu_table(Phy_merged_zhan_raw),tax_table(Phy_merged_zhan_raw),sample_data(zhan_meta_all))





###Correlate PC1 and Shannon

#compute Shannon & Chao index
plot_data_coi=coi_metazoo_meta_all %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID_short") %>%
  mutate(estimate_richness(Phy_merged_coi_raw, measures="Shannon")) %>% 
  mutate(estimate_richness(Phy_merged_coi_raw, measures="Chao1"))

plot_data_18s=zhan_meta_all %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID_short") %>%
  mutate(estimate_richness(Phy_merged_zhan_raw, measures="Shannon")) %>% 
  mutate(estimate_richness(Phy_merged_zhan_raw, measures="Chao1"))



#Correlate
#Test corr
# Calculate Pearson's correlation coefficient and p-value
correlation_result_coi <- cor.test(plot_data_coi$PC1, plot_data_coi$Shannon, method = "spearman")
correlation_result_18s <- cor.test(plot_data_18s$PC1, plot_data_18s$Shannon, method = "pearson")

correlation_result_coi
correlation_result_18s
# Extract Pearson's r and p-value
pearsons_r <- correlation_result_coi$estimate
p_value <- correlation_result_coi$p.value

pearsons_r <- correlation_result_18s$estimate
p_value <- correlation_result_18s$p.value

# Print Pearson's r and p-value
print(paste("Pearson's r:", round(pearsons_r, 3)))
print(paste("p-value:", format(p_value, scientific = FALSE)))


# Run linear regression
lm_model <- lm(Shannon ~ PC1, data = plot_data_18s)

# Summary of linear regression
lm_summary <- summary(lm_model)

# Extracting R-squared and p-value
r_squared <- lm_summary
p_value <- lm_summary$coefficients[2, 4]



#Correlation calculations
lm_model_coi <- lm(Shannon ~ PC1, data = plot_data_coi)

# Summary of linear regression
summary(lm_model_coi)

## Create a scatter plot with regression line, confidence intervals, and color by 'cycle'
#My colors for Cycles
txt_sz=24
my_palette=custom_pallete_all()
coi_plot=ggplot(plot_data_coi, aes(x = PC1, y = Shannon)) +
  geom_point(size=8, aes(, shape = day_night, fill = Sample_ID_short))+
  scale_shape_manual(values = c("day" = 21, "night"=23), name= "Day/Night") +
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 2)) +  # Set fill legend to 2 columns
  theme(legend.title = element_text(face = "bold"))+  # Make sure legend titles are styled consistently  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  coord_cartesian(ylim = c(1.5, 4), xlim = c(min(plot_data_18s$PC1), 7))+
  scale_x_continuous(breaks = seq(-6, 7, by = 2))+
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))+
  # scale_y_continuous(breaks = seq(0, 4, length.out = 10))+
  scale_fill_manual(values = my_palette, name="Station ID")+
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="COI") +
  stat_cor(method="pearson", label.x = 0, label.y = 3.5, size=8)+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 4, label.y = 3.7)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))
coi_plot
saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_all.pdf"), 
    plot = coi_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_all.png"), 
    plot = coi_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

#18s all
my_palette=custom_pallete_all()
zhan_plot=ggplot(plot_data_18s, aes(x = PC1, y = Shannon)) +
  geom_point(size=8, aes(shape = day_night, fill = Sample_ID_short))+
  scale_shape_manual(values = c("day" = 21, "night"=23), name= "Day/Night") +
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 2)) +  # Set fill legend to 2 columns
  theme(legend.title = element_text(face = "bold"))+  # Make sure legend titles are styled consistently  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  coord_cartesian(ylim = c(0, 4), xlim = c(min(plot_data_18s$PC1), 7))+
  scale_x_continuous(breaks = seq(-6, 7, by = 2))+
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))+
  # scale_y_continuous(breaks = seq(0, 4, length.out = 10))+
  scale_fill_manual(values = my_palette, name="Station ID")+
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="18S",
       shape="Cycle", fill="Cycle") +
  stat_cor(method="pearson", label.x = 0, label.y = 3.5, size=8)+
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 4, label.y = 3.7)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title = element_text(size = txt_sz),
        strip.text = element_text(size = txt_sz))
zhan_plot
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_all.pdf"), 
    plot = zhan_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_all.png"), 
    plot = zhan_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  ) }





#Both
both_plot=grid.arrange(coi_plot,zhan_plot)
both_plot

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all.pdf"), 
    plot = both_plot,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all.png"), 
    plot = both_plot,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }


#Save
#COI
output_path <- here("plots", "Q1_physical_analysis")
ggsave(file.path(output_path, "pc1_vs_shannon_coi.png"), coi_plot, width = 10, height = 6, units = "in")

#18s
output_path <- here("plots", "Q1_physical_analysis")
ggsave(file.path(output_path, "pc1_vs_shannon_18s.png"), zhan_plot, width = 10, height = 6, units = "in")

#Both plots 
output_path <- here("plots", "Q1_physical_analysis")
ggsave(file.path(output_path, "pc1_vs_shannon_both.png"), both_plot, width = 10, height = 6, units = "in")






#### PC1 vs Shannon for each size
#compute Shannon index
shannon_coi=estimate_richness(Phy_coi_raw, measures="Shannon") %>% 
  rownames_to_column("Sample_ID")
shannon_18s=estimate_richness(Phy_zhan_raw, measures="Shannon") %>% 
  rownames_to_column("Sample_ID")


plot_data_coi=coi_metazoo_meta %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(.,shannon_coi, by="Sample_ID")%>%
  mutate(max_size=as.factor(max_size))

plot_data_18s=zhan_meta %>% as.data.frame() %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(.,shannon_18s, by="Sample_ID")%>%
  mutate(max_size=as.factor(max_size))

#Test corr
# Calculate Pearson's correlation coefficient and p-value
correlation_result <- cor.test(plot_data_coi$PC1, plot_data_coi$Shannon, method = "pearson")

# Extract Pearson's r and p-value
pearsons_r <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print Pearson's r and p-value
print(paste("Pearson's r:", round(pearsons_r, 3)))
print(paste("p-value:", format(p_value, scientific = TRUE)))

#
facet_correlation <- plot_data_ %>%
  group_by(max_size) %>%
  summarise(pearsons_r = cor(PC1, Shannon, method = "pearson"), 
            p_value = cor.test(PC1, Shannon, method = "pearson")$p.value)

# Print the result
print(facet_correlation)


txt_sz=24
# Mutate the day_night_0_1 column to have "day" and "night" instead of 0 and 1
plot_data_coi <- plot_data_coi %>%
  mutate(day_night = ifelse(day_night_0_1 == 0, "day", "night"))

# Updated ggplot code
coi_plot_sized = ggplot(plot_data_coi, aes(x = PC1, y = Shannon)) +
  geom_point(size=8, aes(, shape = day_night, fill = Sample_ID_short))+
  scale_shape_manual(values = c("day" = 21, "night"=23), name= "Day/Night") +
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 2)) +  # Set fill legend to 2 columns
  scale_fill_manual(values = my_palette) +  # Ensure colors are assigned from the palette
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, aes(color = max_size), show.legend = FALSE, alpha = 0.5) +  # Do not show the color legend for the line
  coord_cartesian(ylim = c(0, 5), xlim = c(min(plot_data_coi$PC1), max(plot_data_coi$PC1))) +
  scale_x_continuous(breaks = seq(-6, max(plot_data_coi$PC1), by = 1.5)) +
  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
  labs(x = "Offshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title = "COI") +
  stat_cor(method = "pearson", label.x = 2, label.y = 4, size = 8) +
  theme_classic() +
  facet_wrap(~max_size, nrow = 3, labeller = labeller(max_size = c("0.5" = "0.2-0.5 mm", "1" = "0.5-1 mm", "2" = "1-2 mm"))) +
  theme(strip.text = element_text(size = txt_sz),
        axis.text.x = element_text(size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title.x = element_text(size = txt_sz),
        axis.title.y = element_text(size = txt_sz)) +
  guides(shape = guide_legend(title = "Day/Night"), fill = guide_legend(title = "Cycle"))  # Only show the legend for points
coi_plot_sized


saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_sized.pdf"), 
    plot = coi_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_coi_sized.png"), 
    plot = coi_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }


#Test corr
# Calculate Pearson's correlation coefficient and p-value
correlation_result <- cor.test(plot_data_18s$PC1, plot_data_18s$Shannon, method = "pearson")

# Extract Pearson's r and p-value
pearsons_r <- correlation_result$estimate
p_value <- correlation_result$p.value

# Print Pearson's r and p-value
print(paste("Pearson's r:", round(pearsons_r, 3)))
print(paste("p-value:", format(p_value, scientific = TRUE)))

#
facet_correlation <- plot_data_18s %>%
  group_by(max_size) %>%
  summarise(pearsons_r = cor(PC1, Shannon, method = "pearson"), 
            p_value = cor.test(PC1, Shannon, method = "pearson")$p.value)

# Print the result
print(facet_correlation)


#Plot to visualize
zhan_plot_sized <- ggplot(plot_data_18s, aes(x = PC1, y = Shannon)) +
  geom_point(size=8, aes(, shape = day_night, fill = Sample_ID_short))+
  scale_shape_manual(values = c("day" = 21, "night"=23), name= "Day/Night") +
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 2)) +  # Set fill legend to 2 columns
  scale_fill_manual(values = my_palette) +  # Ensure colors are assigned from the palette
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, aes(color = max_size), show.legend = FALSE, alpha = 0.5) +  # Do not show the color legend for the line
  coord_cartesian(ylim = c(0, 5), xlim = c(min(plot_data_coi$PC1), max(plot_data_coi$PC1))) +
  scale_x_continuous(breaks = seq(-6, max(plot_data_coi$PC1), by = 1.5)) +
  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
  labs(x = "Offshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title = "COI") +
  stat_cor(method = "pearson", label.x = 2, label.y = 4, size = 8) +
  theme_classic() +
  facet_wrap(~max_size, nrow = 3, labeller = labeller(max_size = c("0.5" = "0.2-0.5 mm", "1" = "0.5-1 mm", "2" = "1-2 mm"))) +
  theme(strip.text = element_text(size = txt_sz),
        axis.text.x = element_text(size = txt_sz),
        axis.text.y = element_text(size = txt_sz),
        axis.title.x = element_text(size = txt_sz),
        axis.title.y = element_text(size = txt_sz)) +
  guides(shape = guide_legend(title = "Day/Night"), fill = guide_legend(title = "Cycle"))  # Only show the legend for points
zhan_plot_sized

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_sized.pdf"), 
    plot = zhan_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_18s_sized.png"), 
    plot = zhan_plot_sized,
    width = 12,  # Width in inches
    height = 10  # Height in inches
  ) }


both_plot_sized=grid.arrange(coi_plot_sized,zhan_plot_sized, nrow=2)
both_plot_sized

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all_sized.pdf"), 
    plot = both_plot_sized,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }

if (saving==1) {
  ggsave(
    filename = here("plots/Q2_diversity_indices/","pc1_vs_shannon_both_all_sized.png"), 
    plot = both_plot,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  ) }





# Q3: Community Composition Using Treemaps: COI, 18S for majority, minority taxa, onshore and offshore ---------------------------------------------------------------------



# COI ---------------------------------------------------------------------

coi_otu=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
coi_taxa=read.csv(here("data/taxa_files/blast_metazoo_coi.csv")) %>% column_to_rownames("Hash")

#Need metadata with offshore_cluster
meta_treemap=read.csv(here(project_path,"data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,clust_group,cycle, max_size)) %>%
  sample_data(.)


OTU = otu_table(as.matrix(coi_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(meta_treemap)
Phy_coi_treemap=phyloseq(OTU,TAX,meta)

#### Transform to Long
phy_merged_long_coi=phyloseq_transform_to_long((Phy_coi_treemap)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species") 

#### Transform to Long
phy_coi_majority=phyloseq_transform_to_long((Phy_coi_treemap)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", NA, Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family)) %>%
  filter((Order %in% c("Calanoida","Euphausiacea")))

phy_coi_minority=phyloseq_transform_to_long((Phy_coi_treemap)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", NA, Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family)) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea")))


## ALL Majority taxa
p_coi_maj=phyloseq_long_treemap_top(phy_coi_majority, Genus,Family ,"COI All",20,colors=NULL, label_group1 = TRUE)
p_coi_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/coi_majority.png"),
  plot = p_coi_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/coi_majority.pdf"),
  plot = p_coi_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


## ALL Minority taxa
p_coi_min=phyloseq_long_treemap_top(phy_coi_minority, Family, Genus ,"COI All",20,colors=NULL, label_group1 = TRUE)
p_coi_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/coi_minority.png"),
  plot = p_coi_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/coi_majority.pdf"),
  plot = p_coi_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)





###By onshore offshore minority coi
off_on=unique(phy_merged_long_coi$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_coi_minority[phy_coi_minority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Species,Genus,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_coi_on_min=list.plots[[1]]
p_coi_on_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_coi_minority.png"),
  plot = p_coi_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_coi_minority.pdf"),
  plot = p_coi_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_coi_off_min=list.plots[[2]]
p_coi_off_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_coi_minority.png"),
  plot = p_coi_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_coi_minority.pdf"),
  plot = p_coi_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


###By onshore offshore majority coi
off_on=unique(phy_merged_long_coi$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_coi_majority[phy_coi_majority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Species,Genus,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_coi_on_maj=list.plots[[1]]
p_coi_on_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_coi_majority.png"),
  plot = p_coi_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_coi_majority.pdf"),
  plot = p_coi_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_coi_off_maj=list.plots[[2]]
p_coi_off_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_coi_majority.png"),
  plot = p_coi_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_coi_majority.pdf"),
  plot = p_coi_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)





# 18s ---------------------------------------------------------------------

#Need metadata with offshore_cluster
OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(meta_treemap)
Phy_18s_treemap=phyloseq(OTU,TAX,meta)

#### Transform to Long
phy_18s_majority=phyloseq_transform_to_long((Phy_18s_treemap)) %>%
  mutate(Genus = ifelse(is.na(Genus), "unknown Genus", Genus)) %>%
  mutate(Genus = ifelse(Genus == "", "unknown Genus", Genus))%>%
  mutate(Family = ifelse(is.na(Family), "unknown Family", Family)) %>%
  filter((Order %in% c("Calanoida","Euphausiacea"))) 


phy_18s_minority=phyloseq_transform_to_long((Phy_18s_treemap)) %>%
  filter(Genus != "Genus")%>%
  filter(Species != "Species")%>%
  mutate(Species = ifelse(Species == "", "NA", Species))%>%
  mutate(Genus = ifelse(Genus == "", NA, Genus)) %>%
  mutate(Genus = ifelse(is.na(Genus), "unknown Genus", Genus)) %>%
  mutate(Family = ifelse(Family == "", NA, Family)) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea")))


## ALL Majority taxa
p_18s_maj=phyloseq_long_treemap_top(phy_18s_majority, Genus,Family ,"18s All",20,colors=NULL, label_group1 = TRUE)
p_18s_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/18s_majority.png"),
  plot = p_18s_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/18s_majority.pdf"),
  plot = p_18s_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


## ALL Minority taxa
p_18s_min=phyloseq_long_treemap_top(phy_18s_minority, Genus,Family ,"18s All",20,colors=NULL, label_group1 = TRUE)
p_18s_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/18s_minority.png"),
  plot = p_18s_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/18s_majority.pdf"),
  plot = p_18s_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)







###By onshore offshore minortiy
off_on=unique(phy_18s_majority$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_18s_minority[phy_18s_minority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Genus,Family,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_18s_on_min=list.plots[[1]]
p_18s_on_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_18s_minority_fam.png"),
  plot = p_18s_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_18s_minority_fam.pdf"),
  plot = p_18s_on_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_18s_off_min=list.plots[[2]]
p_18s_off_min

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_18s_minority_fam.png"),
  plot = p_18s_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_18s_minority_fam.pdf"),
  plot = p_18s_off_min,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


###By onshore offshore majority
off_on=unique(phy_18s_majority$offshore_onshore)
list.plots <- vector('list', length(off_on))
titles=c("Onshore","Offshore")

for (i in 1:length(off_on)){
  phy_sel=phy_18s_majority[phy_18s_majority$offshore_onshore==off_on[i],]
  list.plots[[i]]=phyloseq_long_treemap_top(phy_sel,Species,Genus,titles[i],colors=NULL,top=10, label_group1 = TRUE)
  rm(phy_sel)
}

#Onshore
p_18s_on_maj=list.plots[[1]]
p_18s_on_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/onshore_18s_majority_fam.png"),
  plot = p_18s_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/onshore_18s_majority_fam.pdf"),
  plot = p_18s_on_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)


#Offshore
p_18s_off_maj=list.plots[[2]]
p_18s_off_maj

#PNG & PDF Save
ggsave(
  filename = here("plots/treemaps/offshore_18s_majority_fam.png"),
  plot = p_18s_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/offshore_18s_majority_fam.pdf"),
  plot = p_18s_off_maj,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)




# Using Fido taxa ---------------------------------------------------------------------

# COI ---------------------------------------------------------------------

#Taxa
coi_taxa=read.csv(here("data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv")) %>%
  mutate(Genus = ifelse(Genus == "Genus", Family, Genus)) %>%
  column_to_rownames("Genus") %>% 
  mutate(Hash=X) %>%
  select(-X)



#Predicted proportions
fido_s1_raw=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_taxa_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Genus", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Genus", all = TRUE)%>%
  column_to_rownames("Genus") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_coi_merged_raw


#Make phyloseq objects

#coi
# OTU = otu_table(as.matrix(coi_otu), taxa_are_rows = TRUE)
OTU = otu_table(as.matrix(fido_coi_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(coi_meta)

#USe proportions
phy_coi=phyloseq_transform_to_long((phyloseq(OTU, TAX, meta))) %>% 
  mutate(Genus=asv_code)
phy_coi %>%
  filter(max_size==0.5) %>% 
  phyloseq_long_treemap_top(., Genus, Family ,"",10,colors=NULL, label_group1 = TRUE)->s1_coi
phy_coi %>%
  filter(max_size==1) %>% 
  phyloseq_long_treemap_top(., Genus, Family ,"",10,colors=NULL, label_group1 = TRUE)->s2_coi
phy_coi %>%
  filter(max_size==2) %>% 
  phyloseq_long_treemap_top(., Genus, Family ,"",10,colors=NULL, label_group1 = TRUE)->s3_coi

fido_coi_tree=grid.arrange(s1_coi,s2_coi,s3_coi,nrow=1)

ggsave(
  filename = here("plots/treemaps/fido_coi_treemap.png"),
  plot = fido_coi_tree,
  width = 20,  # Width in inches
  height = 9  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/fido_coi_treemap.pdf"),
  plot = fido_coi_tree,
  width = 20,  # Width in inches
  height = 9  # Height in inches
)



# 18S ---------------------------------------------------------------------



#Predicted proportions
fido_s1_raw=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Family", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Family", all = TRUE)%>%
  column_to_rownames("Family") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_raw

#Metadata
env_metadata_phy=zhan_meta 

#Make phyloseq objects

#Load taxa
zhan_taxa=read.csv(here("data/phyloseq_bio_data/18s/fido_18s_family_tax_table.csv")) %>% 
  column_to_rownames("Family") %>% 
  select(-X)
OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(env_metadata_phy)

#USe proportions
phy_18s=phyloseq_transform_to_long((phyloseq(OTU, TAX, meta))) %>% 
  mutate(Family=asv_code)
phy_18s %>%
  filter(max_size==0.5) %>% 
  phyloseq_long_treemap_top(., Family, Class ,"",10,colors=NULL, label_group1 = TRUE)->s1_18s
phy_18s %>%
  filter(max_size==1.0) %>% 
  phyloseq_long_treemap_top(., Family, Class ,"",10,colors=NULL, label_group1 = TRUE)->s2_18s
phy_18s %>%
  filter(max_size==2) %>% 
  phyloseq_long_treemap_top(., Family, Class ,"",10,colors=NULL, label_group1 = TRUE)->s3_18s

fido_18s_tree=grid.arrange(s1_18s,s2_18s,s3_18s,nrow=1)

ggsave(
  filename = here("plots/treemaps/fido_18s_treemap.png"),
  plot = fido_18s_tree,
  width = 20,  # Width in inches
  height = 9  # Height in inches
)

ggsave(
  filename = here("plots/treemaps/fido_18s_treemap.pdf"),
  plot = fido_18s_tree,
  width = 20,  # Width in inches
  height = 9  # Height in inches
)




# Part II: Stacked Bar Plots -------------------------------------------------------



# COI ---------------------------------------------------------------------

Phy_glom_coi <- Phy_merged_coi %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop))

Phy_glom_coi_minority <- Phy_merged_coi %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop)) %>% 
  ungroup()





#Reorder Calalnoida and Euphausiacea for plotting
Phy_glom_coi %>%
  ungroup() %>% 
  select(Order) %>% 
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>% 
  unique() %>% 
  as.matrix()->other_orders 

Phy_glom_coi=Phy_glom_coi%>%
  mutate(Order=factor(Order, levels =c("Calanoida","Euphausiacea",other_orders)))

#PC1 Labels
labels_for_PC1=Phy_glom_coi %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

#Taxa color maps
orders_coi=unique(Phy_glom_coi$Order)

# Define two Brewer palettes
palette1 <- brewer.pal(8, "Accent")
palette2 <- brewer.pal(8, "Pastel1")

# Join the two palettes
contrast_palette <-  c("#1f77b4", "#ff7f0e", "#2ca02c", "#9edae5" , "#9467bd", "#8c564b", 
                       "#17becf", "#7f7f7f", "#bcbd22", "#e377c2", "#aec7e8", "#ffbb78",
                       "#98df8a", "#ff9896", "#ffaec9","#fc798a")
palette_named <- setNames(contrast_palette, orders_coi)


#Plot all
Phy_glom_coi %>%
  group_by(Order) %>%
  ggplot(.,aes(x = as.factor(PC1), y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Proportion of Total Reads", fill = "Order") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_named) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)-> all_p 



all_p

#Minority plot
Phy_glom_coi %>%
  group_by(Order) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  ggplot(.,aes(x = as.factor(PC1), y = asin(sqrt((prop))), fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Arcsine-Squareroot\nProportion of Residual Reads", fill = "Order") +
  # ggtitle("Raw Relative Read Abundances By Cycle and Family") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_named) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) -> minority_p 
minority_p
stacked_bar_coi=gridExtra::grid.arrange(all_p,minority_p,ncol=1)


#PNG & PDF Save
ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_coi.png"),
  plot = stacked_bar_coi,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)

ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_coi.pdf"),
  plot = stacked_bar_coi,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)




# 18S ---------------------------------------------------------------------
zhan_otu=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))
zhan_taxa=read.csv(here("data/taxa_files/blast_metazoo_18s.csv")) %>% column_to_rownames("Hash")

OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(meta_treemap)
Phy_18s_treemap=phyloseq(OTU,TAX,meta)

Phy_glom_18s <- Phy_18s_treemap %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop))

Phy_glom_18s_minority <- Phy_18s_treemap %>%
  tax_glom(taxrank="Order") %>%
  phyloseq_transform_to_long(.) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  group_by(as.factor(PC1)) %>%
  mutate(prop = n_reads / sum(n_reads),
         total_reads=sum(prop))



# Match categories and get colors for overlapping categories
orders_coi=(unique(Phy_glom_coi$Order))
orders_18s=(unique(Phy_glom_18s$Order))

matching_categories <- intersect(orders_coi,orders_18s)
matching_colors <- palette_named [match(matching_categories, orders_coi)]

# New categories in the second dataset
new_categories <- setdiff(orders_18s, matching_categories)


# Generate new colors for new categories
new_colors <- brewer.pal(length(new_categories), "Pastel2")

# Combine colors for both datasets
palette_18s <- c(setNames(matching_colors, matching_categories), setNames(new_colors, new_categories))



#Reorder Calalnoida and Euphausiacea for plotting
Phy_glom_18s %>%
  ungroup() %>% 
  select(Order) %>% 
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>% 
  unique() %>% 
  as.matrix()->other_orders 

Phy_glom_18s=Phy_glom_18s%>%
  mutate(Order=factor(Order, levels =c("Calanoida","Euphausiacea",other_orders)))


#PC1 Labels
labels_for_PC1=Phy_glom_18s %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))




Phy_glom_18s %>%
  group_by(Order) %>%
  ggplot(.,aes(x = as.factor(PC1), y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Arcsine-Squareroot\nProportion of Residual Reads", fill = "Order") +
  # ggtitle("Raw Relative Read Abundances By Cycle and Family") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_18s) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)-> all_p 


all_p

Phy_glom_18s %>%
  group_by(Order) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>%
  ggplot(.,aes(x = as.factor(PC1), y = asin(sqrt(prop)), fill = Order)) +
  geom_bar(stat = "identity") +
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = "Proportion of Residual Reads", fill = "Order") +
  theme_minimal() +  # Set axis labels
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),  # Increase x-axis label size
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  scale_fill_manual(values = palette_18s) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) -> minority_p 
minority_p
stacked_bar_18s=gridExtra::grid.arrange(all_p,minority_p,ncol=1)


#PNG & PDF Save
ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_18s.png"),
  plot = stacked_bar_18s,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)

ggsave(
  filename = here("plots/Q3_community_composition/stacked_bar_CC_18s.pdf"),
  plot = stacked_bar_18s,
  width = 15,  # Width in inches
  height = 10  # Height in inches
)



### HEAT MAPS
# Aggregating data by Order and PC1
agg_data <- Phy_glom_18s %>%
  group_by(Order) %>%
  filter(!(Order %in% c("Calanoida","Euphausiacea"))) %>% 
  group_by(Order, PC1,Sample_ID_short) %>%
  summarise(n_reads = sum(prop)) %>%
  ungroup()

order_abundance <- agg_data %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(n_reads)) %>%
  arrange((total_abundance)) %>%
  pull(Order)

# Reorder the Order factor based on total abundance
agg_data$Order <- factor(agg_data$Order, levels = order_abundance)

# Creating the heatmap using ggplot
ggplot(agg_data, aes(x = as.factor(PC1), y = Order, fill = asin(sqrt((n_reads))))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +  # Adjust color gradient as needed
  theme_minimal() +
  labs(x = "PC1", y = "Order", fill = "n_reads")+
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability













