#Script to loop through different environmental variables and compare against taxa from quantitative analysis



# Load libraries ----------------------------------------------------------
library(phyloseq)
library(tidyverse)
library(fido)
library(here)



# Raw Reads ---------------------------------------------------------------


# 18S ---------------------------------------------------------------------

#coi read in different sizes
otu18s1=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv")) %>%
  select(starts_with("C"),-X,Family)%>%
  pivot_longer(-Family, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Family) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Family") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otu18s2=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv"))%>%
  select(starts_with("C"),-X,Family)%>%
  pivot_longer(-Family, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Family) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Family") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otu18s3=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv"))%>%
  select(starts_with("C"),-X,Family)%>%
  pivot_longer(-Family, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Family) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Family") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)


zhan_taxa=read.csv(here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv")) %>% 
  select(-X) %>% 
  distinct() %>% 
  mutate(Family2=Family) %>%
  column_to_rownames("Family2") 
zhan_taxa=tax_table(as.matrix(zhan_taxa))



meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(c(-Sizefractionmm,-offshore_onshore,-clust_group,-cycle, -max_size)) %>%
  sample_data(.)

#Taking the mean of replicates
dat_1=phyloseq(otu18s1,zhan_taxa,meta18s) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otu18s2,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otu18s3,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)




set.seed(899)

#Loop thru fido models for each size
counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-mixedlayerdepths,
             -beam_depth,-chl_max,-intergrated_chl,-PAR_1_depth_adj,-day_night_0_1,
             -PC1)) %>% 
    mutate(distance_from_shore=log(distance_from_shore),
           chl_max_depth=log(chl_max_depth),
           nitracline_depth=log(nitracline_depth),
           hypoxia_depth=log(hypoxia_depth))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(formula_obj, data=sample_dat))
  Y <- otu_table(dat) %>% t(.)
  
  #Determine Gamma
  gamma <- c(1,2,3,5,8,10,15,20,50,100,500,700,1000)
  logML <- rep(NA, length(gamma))
  for(i in 1:length(gamma)){
    fit <- pibble(Y, X, Gamma = gamma[i]*diag(nrow(X)), n_samples=5000)
    logML[i] <- fit$logMarginalLikelihood
    print(i)
  }
  plot(gamma, logML, type = "l")
  points(gamma, logML)
  gamma=20
  
  
  ## This is all prior specification
  upsilon <- ntaxa(dat)+3 
  Omega <- diag(ntaxa(dat))
  G <- cbind(diag(ntaxa(dat)-1), -1)
  Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)
  Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
  Gamma <- gamma*diag(nrow(X))
  
  ##This code is used to check priors, not for actual model fitting.
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
  print(priors)
  
  priors <- to_clr(priors)  
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  
  names_covariates(priors) <- rownames(X)
  priors$Y <- Y # remember pibblefit objects are just lists
  posterior <- refit(priors, optim_method="lbfgs", jitter = 1e-5)
  
  tax <- tax_table(dat)[,c("Family")] %>% as.data.frame() %>%
    rownames_to_column("Family2")%>% select(Family)
  num <- 1:nrow(tax)
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,sep="_")
  names_categories(posterior) <- tax
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
  
  #Add size column
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  posterior_summary=posterior_summary %>% 
    mutate(size=sizes[counter])
  assign(paste0("posterior_summary_S", counter), posterior_summary)
  
  ##Let's examine this more
  head(posterior_summary)
  ##Mean is the mean of the posterior samples. You can think of it as the estimated beta for the regression model of that specific taxa.
  ##Covariate: we are fitting y = \beta_0 + \beta_1 * potemp2. So we have estimates both for the intercept and slope.
  ## p2.5, p25, etc. these are the 2.5th, 25th, etc. quantiles of the posterior distribution
  ## NOte that p2.5 and p97.5 would give a 95% credible interval.
  ##So for the first taxa, the intercept has a 95% interval of -2.78,12.5
  ##We assess significance by seeing if zero is in this interval. So, for above, zero is in the inteval, this intercept term isn't significant.
  
  ##Now, we are filtering the posterior summary to significant samples only.
  focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
  focus ##Note that there are zero rows --> no evidence of an effect of potemp2 on any of the taxa
  
  #Plot model wiht fido plot fxn
  fido.p=plot(posterior, par="Lambda", focus.cov = rownames(X)[2:nrow(X)])
  assign(paste0("fido.p_S", counter), fido.p)
  if (nrow(focus) == 0) {
    # Start next loop
    # Your code for the next loop goes here
  } else {
    # Continue with the rest of your code using the 'focus' dataframe
    ##This code will only work if there is a sig. result returned
    focus <- unique(focus$coord)
    # focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p25) == sign(posterior_summary$p75)]]
    focus_cov=rownames(X)
    focus_cov[2]
    
    
    
    
  }
}

grid.arrange(fido.p_S1,fido.p_S2,fido.p_S3, nrow=3)
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)
plot(posterior, par="Lambda", focus.cov = rownames(X)[2:8])

posterior_summary_all=rbind(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3) 


posterior_summary_all%>%
  # filter(sign(p25) == sign(p75))%>%
  filter(covariate != "(Intercept)") %>% 
  # filter(covariate == "PC1") %>% 
  ggplot(., aes(x = mean, y = coord, color=size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = p25, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape=size),size = 3, alpha=0.8) + 
  facet_wrap(~covariate, nrow=4)+
  labs(title="Pibble Model for 18S Taxa",x = "Centered Log-Ratio",
       y="", color = "Size Fraction (mm)", shape="Size Fraction (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16),# Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14))->all_plot_18s_rra
all_plot_18s_rra

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_18s_taxa.pdf"), 
    plot = all_plot_18s_rra,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_18s_taxa.png"), 
    plot = all_plot_18s_rra,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  )}




# COI ---------------------------------------------------------------------
#Using agglomerated taxa from fido
otucoi1=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_taxa_phy.csv")) %>%
  select(starts_with("C"),-X,Genus)%>%
  pivot_longer(-Genus, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Genus) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Genus") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otucoi2=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv"))%>%
  select(starts_with("C"),-X,Genus)%>%
  pivot_longer(-Genus, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Genus) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Genus") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otucoi3=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv"))%>%
  select(starts_with("C"),-X,Genus)%>%
  pivot_longer(-Genus, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Genus) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Genus") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)


coi_taxa=read.csv(here("data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv")) %>%
  mutate(Genus = ifelse(Genus == "Genus", Family, Genus)) %>%
  column_to_rownames("Genus") %>% 
  mutate(Genus=row.names(.)) %>%
  mutate(Hash=X) %>%
  select(-X,-Hash) %>%
  select(1:8, 10, 9) 
coi_taxa=tax_table(as.matrix(coi_taxa))



metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,cycle, max_size)) %>%
  sample_data(.)

dat_1=phyloseq(otucoi1,coi_taxa,metacoi) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otucoi2,coi_taxa,metacoi)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otucoi3,coi_taxa,metacoi)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)

set.seed(899)


### =========== PC1
data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)


#Check distribution of vars
# Convert the data to long format
metacoi_long <- metacoi %>% 
  gather(key = "Variable", value = "Value") %>%
  distinct() %>% 
  filter(!is.na(as.numeric(Value)))
  

# Create histograms and facet them
ggplot(metacoi, aes(x = PC1)) +   
  geom_histogram(bins=5, color = "black", fill = "lightblue", alpha = 0.6) +  
  # facet_wrap(~ Variable, scales = "free") +
  labs(title = "Histograms of All Columns", x = "Values", y = "Frequency")

set.seed(899)


counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-mixedlayerdepths,
             -beam_depth,-chl_max,-intergrated_chl,-PAR_1_depth_adj,-day_night_0_1,
             -PC1)) %>% 
    mutate(distance_from_shore=log(distance_from_shore),
           chl_max_depth=log(chl_max_depth),
           nitracline_depth=log(nitracline_depth),
           chl_max_depth=log(chl_max_depth))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(formula_obj, data=sample_dat))
  Y <- otu_table(dat) %>% t(.)
  
  
  
  ## This is all prior specification
  upsilon <- ntaxa(dat)+3 
  Omega <- diag(ntaxa(dat))
  G <- cbind(diag(ntaxa(dat)-1), -1)
  Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)
  Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
  Gamma <- gamma*diag(nrow(X))
  
  ##This code is used to check priors, not for actual model fitting.
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
  print(priors)
  
  priors <- to_clr(priors)  
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  
  names_covariates(priors) <- rownames(X)
  priors$Y <- Y # remember pibblefit objects are just lists
  posterior <- refit(priors, optim_method="lbfgs", jitter = 1e-5)
  plot(posterior)
  tax <- tax_table(dat)[,c("Genus")] %>% as.data.frame() %>%
    rownames_to_column("Genus2")%>% select(Genus)
  num <- 1:nrow(tax)
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,sep="_")
  names_categories(posterior) <- tax
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
  
  #Add size column
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  posterior_summary=posterior_summary %>% 
    mutate(size=sizes[counter])
  assign(paste0("posterior_summary_S", counter), posterior_summary)
  
  ##Let's examine this more
  head(posterior_summary)
  ##Mean is the mean of the posterior samples. You can think of it as the estimated beta for the regression model of that specific taxa.
  ##Covariate: we are fitting y = \beta_0 + \beta_1 * potemp2. So we have estimates both for the intercept and slope.
  ## p2.5, p25, etc. these are the 2.5th, 25th, etc. quantiles of the posterior distribution
  ## NOte that p2.5 and p97.5 would give a 95% credible interval.
  ##So for the first taxa, the intercept has a 95% interval of -2.78,12.5
  ##We assess significance by seeing if zero is in this interval. So, for above, zero is in the inteval, this intercept term isn't significant.
  
  ##Now, we are filtering the posterior summary to significant samples only.
  focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
  focus ##Note that there are zero rows --> no evidence of an effect of potemp2 on any of the taxa
  
  if (nrow(focus) == 0) {
    # Start next loop
    # Your code for the next loop goes here
  } else {
    # Continue with the rest of your code using the 'focus' dataframe
    ##This code will only work if there is a sig. result returned
    focus <- unique(focus$coord)
    focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p25) == sign(posterior_summary$p75)]]
    focus_cov=rownames(X)
    focus_cov[2]
    
    
    
    
  }
}


posterior_summary_all=rbind(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3) 

posterior_summary_all%>% 
  filter(covariate != "(Intercept)") %>% 
  # filter(covariate == "PC1") %>% 
  ggplot(., aes(x = mean, y = coord, color=size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = p25, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape=size),size = 3, alpha=0.8) + 
  facet_wrap(~covariate, nrow=3)+
  labs(title="Pibble Model for COI Taxa",x = "Centered Log-Ratio",
       y="", color = "Size Fraction (mm)", shape="Size Fraction (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) ->all_plot_coi

all_plot_coi

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_coi_taxa_pcr_ra.pdf"), 
    plot = all_plot_coi,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_coi_taxa_pcr_ra.png"), 
    plot = all_plot_coi,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}




# PCR-RA (not working yet) ------------------------------------------------------------------

# 18S ---------------------------------------------------------------------

here()
#coi read in different sizes
otu18s1_pcrra=read.csv(here("data/predicted_og/predicted_og_18s_04_04_2024_s1_phy.csv")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = gsub('predicted ', '', replicate)) %>% 
  select(n_reads,coord,Sample_ID) %>% 
  pivot_wider(names_from = Sample_ID, values_from = n_reads) %>% 
  column_to_rownames("coord") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

# otu18s1_pcrra_counts <- otu_table(otu18s1) * otu_table(otu18s1_pcrra)


otu18s2_pcrra=read.csv(here("data/predicted_og/predicted_og_18s_04_04_2024_s1_phy.csv")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = gsub('predicted ', '', replicate)) %>% 
  select(n_reads,coord,Sample_ID) %>% 
  pivot_wider(names_from = Sample_ID, values_from = n_reads) %>% 
  column_to_rownames("coord") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)
otu18s2_pcrra_counts <- otu_table(otu18s2) * otu_table(otu18s2_pcrra)

otu18s3_pcrra=read.csv(here("data/predicted_og/predicted_og_18s_04_04_2024_s1_phy.csv")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = gsub('predicted ', '', replicate)) %>% 
  select(n_reads,coord,Sample_ID) %>% 
  pivot_wider(names_from = Sample_ID, values_from = n_reads) %>% 
  column_to_rownames("coord") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)
otu18s3_pcrra_counts <- otu_table(otu18s3) * otu_table(otu18s3_pcrra)


zhan_taxa=read.csv(here("data/phyloseq_bio_data/18s/fido_18s_family_tax_table.csv")) %>%
  mutate(Family2=Family) %>%
  column_to_rownames("Family2") %>% 
  select(-X)
zhan_taxa=tax_table(as.matrix(zhan_taxa))



meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(c(-Sizefractionmm,-offshore_onshore,-clust_group,-cycle, -max_size)) %>%
  sample_data(.)

dat_1=phyloseq(otu18s1_pcrra_counts,zhan_taxa,meta18s) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otu18s2_pcrra_counts,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otu18s3_pcrra_counts,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)




set.seed(899)


counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-mixedlayerdepths,
             -beam_depth,-chl_max,-intergrated_chl,-PAR_1_depth_adj,-day_night_0_1,
             PC1))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(formula_obj, data=sample_dat))
  Y <- otu_table(dat) %>% t(.)
  
  
  
  ## This is all prior specification
  upsilon <- ntaxa(dat)+3 
  Omega <- diag(ntaxa(dat))
  G <- cbind(diag(ntaxa(dat)-1), -1)
  Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)
  Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
  Gamma <- diag(nrow(X))
  
  ##This code is used to check priors, not for actual model fitting.
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
  print(priors)
  
  priors <- to_clr(priors)  
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  
  names_covariates(priors) <- rownames(X)
  priors$Y <- Y # remember pibblefit objects are just lists
  posterior <- refit(priors, optim_method="lbfgs", jitter = 1e-5)
  
  tax <- tax_table(dat)[,c("Family")] %>% as.data.frame() %>%
    rownames_to_column("Family2")%>% select(Family)
  num <- 1:nrow(tax)
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,sep="_")
  names_categories(posterior) <- tax
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
  
  #Add size column
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  posterior_summary=posterior_summary %>% 
    mutate(size=sizes[counter])
  assign(paste0("posterior_summary_S", counter), posterior_summary)
  
  ##Let's examine this more
  head(posterior_summary)
  ##Mean is the mean of the posterior samples. You can think of it as the estimated beta for the regression model of that specific taxa.
  ##Covariate: we are fitting y = \beta_0 + \beta_1 * potemp2. So we have estimates both for the intercept and slope.
  ## p2.5, p25, etc. these are the 2.5th, 25th, etc. quantiles of the posterior distribution
  ## NOte that p2.5 and p97.5 would give a 95% credible interval.
  ##So for the first taxa, the intercept has a 95% interval of -2.78,12.5
  ##We assess significance by seeing if zero is in this interval. So, for above, zero is in the inteval, this intercept term isn't significant.
  
  ##Now, we are filtering the posterior summary to significant samples only.
  focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
  focus ##Note that there are zero rows --> no evidence of an effect of potemp2 on any of the taxa
  
  if (nrow(focus) == 0) {
    # Start next loop
    # Your code for the next loop goes here
  } else {
    # Continue with the rest of your code using the 'focus' dataframe
    ##This code will only work if there is a sig. result returned
    focus <- unique(focus$coord)
    # focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p25) == sign(posterior_summary$p75)]]
    focus_cov=rownames(X)
    focus_cov[2]
    
    
    
    
  }
}


posterior_summary_all=rbind(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3) 

posterior_summary_all%>% 
  filter(covariate != "(Intercept)") %>% 
  # filter(covariate == "PC1") %>% 
  ggplot(., aes(x = mean, y = coord, color=size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = p25, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape=size),size = 3, alpha=0.8) + 
  facet_wrap(~covariate, nrow=6)+
  labs(title="Pibble Model for 18S Taxa",x = "Centered Log-Ratio",
       y="", color = "Size Fraction (mm)", shape="Size Fraction (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) ->all_plot_18s_pcrra
all_plot_18s_pcrra

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_18s_taxa_pcr_ra.pdf"), 
    plot = all_plot_18s_pcrra,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_18s_taxa_pcr_ra.png"), 
    plot = all_plot_18s_pcrra,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}


args_null=function (par, argl, default) {
  if (is.null(argl[[par]])) 
    return(default)
  return(argl[[par]])
}
