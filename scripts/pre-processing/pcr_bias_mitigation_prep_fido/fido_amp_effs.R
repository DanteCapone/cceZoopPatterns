library(fido)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(5903)

# 18S ---------------------------------------------------------------------


###Load in the ECDF-filtered data for the 18S primer using long format species and hash name so I ca identify taxa
##First Size 1
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/meta_18s_unaveraged_s1.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 


# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.2-0.5 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
##First Size 1
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/meta_18s_unaveraged_s2.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3
##First Size 1
#Phyloseq Filtered
fido_input_filt_s3=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s_s3=read.csv(file.path("data/fido/meta_18s_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt_s3))
colnames(fido_input_filt_s3) <- gsub("^X", "", colnames(fido_input_filt_s3))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s_s3 <- meta_18s_s3[match(colnames(fido_input_filt_s3), meta_18s_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size

p_s3
#Plot all 
# Combine legends
amp_effs_18s=grid.arrange(p_s1,p_s2,p_s3)


#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s.png"),
  plot = amp_effs_18s,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s.pdf"),
  plot = amp_effs_18s,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)






# COI ---------------------------------------------------------------------
###Load in the ECDF-filtered data for the 18S primer using long format species and hash name so I ca identify taxa
##First Size 1
#Phyloseq Filtered


fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_s1.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 


# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.2-0.5 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
##First Size 1
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_s2.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3
##First Size 1
#Phyloseq Filtered
fido_input_filt_s3=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi_s3=read.csv(file.path("data/fido/meta_coi_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt_s3))
colnames(fido_input_filt_s3) <- gsub("^X", "", colnames(fido_input_filt_s3))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi_s3 <- meta_coi_s3[match(colnames(fido_input_filt_s3), meta_coi_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3

#Plot all 
# Combine legends
amp_effs_coi=grid.arrange(p_s1,p_s2,p_s3)


#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi.png"),
  plot = amp_effs_coi,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi.pdf"),
  plot = amp_effs_coi,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)

