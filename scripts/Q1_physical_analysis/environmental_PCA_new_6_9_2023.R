#1) load in metadata and select desired rows 

library(knitr)
library(tidyverse)
library(clustertend)
library(cluster)
library(flexclust)
library(fpc)
library(clustertend)
library(ClusterR)
library(factoextra)
library(gridExtra)
library(paran)
library(BBmisc)
library(stats)
library(car)
library(ggfortify)
library(ggrepel)
library(here)
library(missMDA)
library(factoextra)
library(NbClust)


#Set working directory
here()
env_metadata<-read.csv(here("data/")) 
names(env_metadata)

#temp
#density
#chl max depth
#chla 
#dist from shore

#sal
#"oxy_sat"
#max beam depth
#chl max value
#hypoxia depth
#MLD
#1% Par depth
#integrated chl
#nitracline
#NO2
#NO3
#PO4
#SIL
#NH4
selected_vars <- c(Sample_ID)
selected_vars <- c(Sample_ID,Sample_ID_short,Sample_ID_dot,potemp2,density2,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                   PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,integrated_chl,day_night_0_1,distance_from_shore) 

env_metadata_long= env_metadata %>% dplyr::select(c(Sample_ID,max_size,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                   PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                   PAR_1_depth)) 
# write.csv(env_metadata_long,"data/CURRENT_WORKING_Metadata/env_metadata_all_samples_8.16.2023.csv")

env_metadata_phy= env_metadata %>% dplyr::select(c(Sample_ID_dot,Sample_ID,max_size,Sizefractionmm,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                    PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                    PAR_1_depth)) 
env_metadata_phy_map= env_metadata %>% dplyr::select(c(Sample_ID_dot,Sample_ID,max_size,Sizefractionmm,Sample_ID_short,cycle,potemp2,density2,oxy_sat,Sal2,NO2,nitracline_depth,NO3,fluorescence,hypoxia_depth,chl_max_depth,
                                                   PO4,SIL,NH4,beam_depth,chl_max,mixedlayerdepths,intergrated_chl,day_night_0_1,distance_from_shore,
                                                   PAR_1_depth,LONGITUDE,LATITUDE)) 
# env_metadata_phy=env_metadata_phy_map

#If I want to just look at sampling site
env_metadata_sel= env_metadata_long %>% filter(max_size==2.0) %>% dplyr::select(-c(max_size))
write.csv(env_metadata_sel,"data/CURRENT_WORKING_Metadata/env_metadata_just_sites_8.16.2023.csv")

#Else
env_metadata_sel= env_metadata_long


#Make nighttime PAR equal to day time
env_metadata_sel =env_metadata_sel %>%
  group_by(cycle) %>%
  mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>%
  filter(max_size==0.5) %>%
  dplyr::select(-max_size)%>% 
  bind_rows() %>% ungroup() %>%
  dplyr::select(-c("Sample_ID","cycle","PAR_1_depth"))%>%column_to_rownames("Sample_ID_short") 




## Imputation
estim_ncpPCA(env_metadata_sel)
metadata_impute <- imputePCA(env_metadata_sel,ncp=1)





#####BIOLOGICAL METADATA Make an imputed dataset for phyloseq analysis

env_metadata_phy_sel =env_metadata_phy %>%
  group_by(cycle) %>%
  mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>% 
  bind_rows() %>% ungroup() %>%
  dplyr::select(-c("Sample_ID_short","cycle","PAR_1_depth","max_size","Sizefractionmm","Sample_ID_dot"))%>%column_to_rownames("Sample_ID")

#For cluster analysis use short DF (8/8/2023)
env_metadata_phy_sel =env_metadata_phy %>%
  group_by(cycle) %>%
  mutate(PAR_1_depth_adj = ifelse(cycle %in% c(1,2,3), max(PAR_1_depth),PAR_1_depth)) %>% 
  bind_rows() %>% ungroup() %>% 
  dplyr::select(-c("Sample_ID","cycle","PAR_1_depth","max_size","Sizefractionmm","Sample_ID_dot"))%>%
  distinct()%>%column_to_rownames("Sample_ID_short")


## Imputation
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

# 8/7/2023 Save
# write.csv(metadata_impute_df_phy,"data/CURRENT_WORKING_Metadata/env_metadata_impute_for_clusters_8.7.2023.csv")


# write.csv(metadata_impute_df_phy,"data/CURRENT_WORKING_Metadata/env_metadata_impute_phyloseq_5.22.2023.csv")
########


#Saving
#write.csv(metadata_impute_df,"data/CURRENT_WORKING_Metadata/env_metadata_impute_5.22.2023.csv")

#PCA
pca_in=metadata_impute$completeObs %>% as.data.frame() 
env_pca=prcomp(pca_in, scale=TRUE)
# env_pca=PCA(pca_in)
summary(env_pca)
env_pca_df=env_pca$x %>% as.data.frame()



#New plot
# Extract cos^2 values for variables on Dim.1 and Dim.2
var_cos2 <- env_pca$var$cos2
ggplot(var_cos2, aes(x = Dim.1)) + 
  geom_histogram(fill = "blue", color = "black", alpha = 0.7, bins = 30) + 
  labs(title = "Histogram of cos^2 values for Dim.1",
       x = "cos^2 values",
       y = "Frequency") +
  theme_minimal()

var_cos2 %>% as.data.frame() %>%
  dplyr::select(Dim.1, Dim.2) %>%
  gather(key = "Dimension", value = "cos2") %>%
  ggplot(aes(x = Dimension, y = cos2)) +
  geom_boxplot(fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Boxplot of cos^2 values for the first two dimensions",
       x = "Dimensions",
       y = "cos^2 values") +
  theme_minimal()

# Determine a threshold for significance. Here, we use a threshold of 0.6 as an example.
threshold <- 4

# Filter out variables that don't meet the threshold on both dimensions
significant_vars <- rownames(var_cos2)[var_cos2[,1] > threshold | var_cos2[,2] > threshold]

# Use fviz_pca_var() to plot only significant variables
# We'll utilize the `select.var` argument to select the significant variables.
fviz_pca_var(env_pca, select.var = list(name = significant_vars), repel = TRUE)

##Contribution of vars to pcs

#Which PCs contribute significantly
explained_var <- env_pca$sdev^2 / sum(env_pca$sdev^2)
df <- data.frame(component = 1:length(explained_var), explained_var = explained_var)
fviz_eig(env_pca, addlabels = TRUE)

var<-get_pca_var(env_pca)
a<-fviz_contrib(env_pca,"var",axes = 1) # default angle=45?
plot(a,main = "Variables percentage contribution of PC1")


var<-get_pca_var(env_pca)
a<-fviz_contrib(env_pca,"var",axes = 2) # default angle=45?
plot(a,main = "Variables percentage contribution of PC2")
grid


#Using ggfortify

#Biplot
fviz_pca_biplot(env_pca, 
                pointsize = 3,  # size of data points
                col.var = "cycle", # color of variable vectors
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

#With loadings
pca_2=autoplot(env_pca, data=metadata_impute_df_phy, colour="cycle",label=TRUE, label.size=4, loadings=TRUE,size=3, loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size =4, scale = 0)+theme_classic()
pca_2


grid.arrange(plot1, pca_2, ncol=2)
  


#Clustering
library(pvclust)
set.seed(123)
env_pca_df=env_pca$x %>% as.data.frame()
# env_pca_df=pca_in
res.pv <- pvclust(t(env_pca_df), method.dist="cor",method.hclust="average", nboot = 10000)
# seplot(res.pv, identify=TRUE)

#Plot
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


#Convert res.pv to dataframe and save

write.csv(t(env_pca_df),"data/clustering_analysis/env_pca_dist_for_clustering_8.8.2023.csv")



res.pv.phys=res.pv %>% as.data.frame(.)

#Optimal cluster #
#Elbow Method
fviz_nbclust(pca_in, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(pca_in, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(pca_in, kmeans, nstart = 25,  method = "gap_stat", nboot = 100)+
  labs(subtitle = "Gap statistic method")


#Nbclust
NbClust(data = pca_in, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")



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

write.csv(metadata_impute_df_phy,"data/CURRENT_WORKING_Metadata/env_metadata_impute_phyloseq_6.9.2023.csv")














######SCRAP

#potemp density scatterplot
ggplot(eDNAleray_metazoo_meta, aes(potemp2,density2)) +
  geom_point()+
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()



#Add rownames
row.names(env_metadata_use)=env_metadata_uniq$Sample_ID


num_cols <- unlist(lapply(env_metadata_use, is.numeric))         # Identify numeric columns

env_pca=prcomp(env_metadata_use, scale=TRUE)
summary(env_pca)


var_explained <- env_pca$sdev^2/sum(env_pca$sdev^2)

#Preprocess
env_pca_df=env_pca$x %>% 
  as.data.frame
row.names(env_pca_df)=env_metadata_uniq$Sample_ID
env_pca_df$Cycle=env_metadata_uniq$cycle
env_pca_df$size_fraction=env_metadata_uniq$Size.fraction..mm.

##Eigenvalues
eig.val<-get_eigenvalue(env_pca)
eig.val

#retain 5 pc's based on Kaiser Criterion

##COntribution of vars to pcs
var<-get_pca_var(env_pca)
a<-fviz_contrib(env_pca,"var",axes = 1) # default angle=45?
plot(a,main = "Variables percentage contribution of PC1")


var<-get_pca_var(env_pca)
a<-fviz_contrib(env_pca,"var",axes = 2) # default angle=45?
plot(a,main = "Variables percentage contribution of PC2")

####Vizualize PCAs
library(ggfortify)


#Using ggfortify
#Without loadings
autoplot(env_pca, data=env_metadata_uniq, colour="cycle",shape="day_night",size=3, labels="Sample_ID")+theme_classic()


#With loadings
plot2=autoplot(env_pca, data=env_metadata_uniq, colour="cycle",label=TRUE, label.size=2, loadings=TRUE,size=3, loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size =3)+theme_classic()
plot2
## Ellipses

plot1=autoplot(env_pca, data=env_metadata_uniq, colour="cycle",label=TRUE,labels="Sample_ID",label.size=3)+theme_classic()
plot1

grid.arrange(plot1, plot2, ncol=2)




#####Clustering
library(pvclust)

set.seed(123)
res.pv <- pvclust(t(env_pca$x), method.dist="cor", 
                  method.hclust="average", nboot = 100)
# Default plot
plot(res.pv, hang = -1, cex = 0.5)
pvrect(res.pv)

clusters <- pvpick(res.pv)
clusters
offshore=array(unlist(clusters$clusters[2]),dim=c(8,1))
onshore=array(unlist(clusters$clusters[1]),dim=c(6,1))



#Add clustering variable to dataframe

env_metadata$clust_group=rep(1,length(env_metadata$Sample_ID))
env_metadata$clust_group[env_metadata$Sample_ID %in% offshore]=2
env_metadata$offshore_onshore[env_metadata$clust_group==2]="offshore"
env_metadata$offshore_onshore[env_metadata$clust_group==1]="onshore"
# env_metadata=env_metadata[,-c(37)]
write.csv(env_metadata,"Env_MetaData_clusters_v2_T2.csv")

#Add same to short metadata vector
env_metadata_uniq$clust_group=rep(1,length(env_metadata_uniq$Sample_ID))
env_metadata_uniq$clust_group[env_metadata_uniq$Sample_ID %in% offshore]=2
env_metadata_uniq$offshore_onshore[env_metadata_uniq$clust_group==2]="offshore"
env_metadata_uniq$offshore_onshore[env_metadata_uniq$clust_group==1]="onshore"


#Add PC1 to the metadata file
pc1=env_pca_df[,1]
env_pca_df$Sample_ID=row.names(env_pca_df)
row.names(pc1)
# env_metadata$pca_1=rep(1,length(env_metadata$Sample_ID))
env_metadata$clust_group[env_metadata$Sample_ID %in% offshore]=2
temp=left_join(env_metadata, env_pca_df, by="Sample_ID")

env_metadata=temp[,1:41]
write.csv(env_metadata,"Env_MetaData_clusters_v2_T2.csv")
