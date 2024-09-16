#Script for EDA visualization, pca and clustering

#1) load in metadata and select desired rows 

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


#Set working directory
here()
env_metadata<-read.csv(here("data/pre_processing/metadata03-28-2023.csv")) 
names(env_metadata)


# Making the metadata -----------------------------------------------------

#Various metadata frames for each analysis
env_metadata_long= env_metadata %>% dplyr::select(c(Sample_ID,max_size,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                    PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                    PAR_1_depth)) 
# env_metadata_phy= env_metadata %>% dplyr::select(c(Sample_ID_dot,Sample_ID,max_size,Sizefractionmm,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
#                                                     PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
#                                                     PAR_1_depth)) 
# env_metadata_phy_map= env_metadata %>% dplyr::select(c(Sample_ID_dot,Sample_ID,max_size,Sizefractionmm,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
#                                                    PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
#                                                    PAR_1_depth,LONGITUDE,LATITUDE)) 

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
# env_metadata_phy_sel =env_metadata_phy %>%
#   group_by(cycle) %>%
#   mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>% 
#   bind_rows() %>% ungroup() %>% 
#   dplyr::select(-c("Sample_ID","cycle","PAR_1_depth","max_size","Sizefractionmm","Sample_ID_dot"))%>%
#   distinct()%>%column_to_rownames("Sample_ID_short")




# IMPUTATION --------------------------------------------------------------

## Imputation, impute the missing nutrient data
estim_ncpPCA(env_metadata_sel)
metadata_impute <- imputePCA(env_metadata_sel,ncp=1)

## Imputation on phyloseq metadata frames
estim_ncpPCA(env_metadata_phy_sel)
metadata_impute <- imputePCA(env_metadata_phy_sel,ncp=1)

env_metadata_phy_add_back =env_metadata_phy %>%
  filter(max_size==0.5)
metadata_impute_df_phy=metadata_impute$completeObs %>% data.frame() %>% 
  #Add back in categoraical vars
  mutate(Sizefractionmm=env_metadata_phy_add_back$Sizefractionmm) %>%
  mutate(max_size=env_metadata_phy_add_back$max_size) %>%
  mutate(cycle=env_metadata_phy_add_back$cycle) %>%
  mutate(Sample_ID_short=env_metadata_phy_add_back$Sample_ID_short) %>%
  mutate(Sample_ID_dot=env_metadata_phy_add_back$Sample_ID_dot)   # column_to_rownames("Sample_ID_short")












# Correlation plot --------------------------------------------------------

meta_data_phy=read.csv(here("Zoop_Patterns/data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size))


#Load in the phyloseq data and format to a table

#COI raw reads
otucoi=read.csv(here("Zoop_Patterns/data/phyloseq_bio_data/COI/metazooprunedcoi_otu.csv")) %>%
  column_to_rownames("Hash") %>%
  select(where(~ !is.na(.[[1]])))

taxcoi=read.csv(here("Zoop_Patterns/data/phyloseq_bio_data/COI/metazooprunedcoi_tax.csv")) %>% column_to_rownames("Hash")


dat_all=phyloseq(otucoi,taxcoi,meta_data_phy)
dat=merge_samples(dat_all,"Sample_ID_short",fun= mean)%>%
  filter_taxa(function(x) sum(x > 3) > 0.10*length(x), TRUE)

set.seed(899)


##Env Correlation matrix
meta_corr=meta_data_phy %>% dplyr::select(-Sample_ID_short) %>%
  cor(.)

# Compute p-values using correlation matrix
p_values <- cor.mtest(meta_data_phy %>% dplyr::select(-Sample_ID_short))$p %>%
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
library(corrplot)
dev.off()
# Sort the row names and column names to ensure alignment
meta_corr <- meta_corr[order(rownames(meta_corr)), order(colnames(meta_corr))]
p_vals_adj <- p_vals_adj[order(rownames(p_vals_adj)), order(colnames(p_vals_adj))]


corr_plot=corrplot(meta_corr,p.mat=p_vals_adj, type = 'lower', order = 'FPC', tl.col = 'black',
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
plot(a,main = "Variables percentage contribution of PC1")


var<-get_pca_var(env_pca)
pc2_p<-fviz_contrib(env_pca,"var",axes = 2) # default angle=45?
plot(a,main = "Variables percentage contribution of PC2")
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
#Without loadings
autoplot(env_pca, data=metadata_impute_df_phy,size=3, labels="Sample_ID_dot")+theme_classic()

## Ellipses

plot1=autoplot(env_pca, data=metadata_impute_df_phy,label=TRUE,colour="cycle", label.size=4, loadings=TRUE,size=3, loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size =5, scale = 1, repel=TRUE)+theme_classic()
plot1
  


#Clustering using pvclust. Cluster the PCA dataframe using hierarchical clustering
# 
library(pvclust)
set.seed(123)
env_pca_df=env_pca$x %>% as.data.frame()
# env_pca_df=pca_in
res.pv <- pvclust(t(env_pca_df), method.dist="cor",method.hclust="average", nboot = 1000)
# seplot(res.pv, identify=TRUE)

#Plot and save figures
here()
pdf(here("figures/physiical_clusters_2_2024/pca_cluster_2_28_2024.pdf"),
    width = 16,  # Width in inches
    height = 12)  # Height in inches
plot(res.pv, hang = -1, cex = 1.2,xlab="Sample ID", main="Environmental PCA Clusterings \n with p-values (%)")
pvrect(res.pv, alpha = 0.95, lwd=2)
dev.off()

#PNG
png(here("figures/physiical_clusters_2_2024/pca_cluster_2_28_2024.png"),
    width = 1600,  # Width in pixels
    height = 1200,  # Height in pixels
    res = 300)  # Resolution in dots per inch (dpi)
plot(res.pv, hang = -1, cex = 1,xlab="Sample ID", main="Environmental PCA Clusterings \n with p-values (%)")
pvrect(res.pv, alpha = 0.95, lwd=2)
dev.off()




res.pv.phys=res.pv %>% as.data.frame(.)

#Determine the Optimal cluster #
# Silhouette method
fviz_nbclust(pca_in, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")


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
env_pca_df =env_pca_df%>% rownames_to_column("Sample_ID_short")
row.names(pc1)

metadata_impute_df_phy=left_join(metadata_impute_df_phy, env_pca_df, by="Sample_ID_short") %>% 
  dplyr::select(-tail(names(.), 16))

# write.csv(metadata_impute_df_phy,"data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")











