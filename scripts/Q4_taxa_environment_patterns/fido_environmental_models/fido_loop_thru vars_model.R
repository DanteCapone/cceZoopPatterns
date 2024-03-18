#Script to loop through different environmental variables and compare against taxa from quantitative analysis



library(MicrobeDS)
library(phyloseq)
library(tidyverse)
library(fido)
library(corrplot)
library(here)

here()
#coi read in different sizes
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


set.seed(899)


counter=0
for (dat_name in names(data_list)) {

  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-nitracline_depth,-mixedlayerdepths,-chl_max_depth,
             -hypoxia_depth,-beam_depth,-chl_max,-intergrated_chl,-distance_from_shore,-PAR_1_depth_adj,-day_night_0_1,-density2))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(~PC1, data=sample_dat))
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
  
  tax <- tax_table(dat)[,c("Genus")] %>% as.data.frame() %>%
    rownames_to_column("Genus2")%>% select(Genus)
  num <- 1:nrow(tax)
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,sep="_")
  names_categories(posterior) <- tax
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
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
    focus_cov
    
    # Create the PDF filename string
    save_path="plots/models_by_var/"
    filename <- paste0(save_path,"S_", dat_name, "_PC1_coi", ".png")
    
    # For demonstration purposes, print the filename
    cat("Saving to:", filename, "\n")
    
    
    png(filename, width = 1000, height = 600)
    
    # Your plotting code
    sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
    pp=plot(posterior, par="Lambda", focus.coord = focus, focus.cov = "PC1")+
      labs(title=paste(sizes[counter], "PC1"))
    base::print(pp)
    Sys.sleep(2)
    # Close the PDF device
    dev.off()
    #Save loop plot
    assign(paste0("pp", counter), pp) 
  }
}
all_plot=grid.arrange(pp1,pp2,pp3, nrow=1,
                      bottom = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore")
all_plot



#Now make into a ggplot
#####
#Add size column
posterior_summary_S1 =posterior_summary_S1 %>% mutate(size="0.2-0.5")
posterior_summary_S2 =posterior_summary_S2 %>% mutate(size="0.5-1")
posterior_summary_S3 =posterior_summary_S3 %>% mutate(size="1-2")

posterior_summary_all=bind_rows(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3)%>%
  filter(covariate != "(Intercept)")

#Plot
all_plot=ggplot(posterior_summary_all, aes(x = mean, y = coord, color=size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = p2.5, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape=size),size = 8, alpha=0.8) + 
  labs(title="Pibble Model for COI Taxa vs. PC1",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Size Fraction (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 




saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q3_taxa_vars/fido_PC1_coi_taxa.pdf"), 
    plot = all_plot,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q3_taxa_vars/fido_PC1_coi_taxa.png"), 
    plot = all_plot,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}




### =========== PC1 18S

here()
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


zhan_taxa=read.csv(here("data/phyloseq_bio_data/18s/fido_18s_family_tax_table.csv")) %>%
  mutate(Family2=Family) %>%
  column_to_rownames("Family2") %>% 
  select(-X)
zhan_taxa=tax_table(as.matrix(zhan_taxa))



meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,cycle, max_size)) %>%
  sample_data(.)

dat_1=phyloseq(otu18s1,zhan_taxa,meta18s) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otu18s2,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otu18s3,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)




set.seed(899)


counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-nitracline_depth,-mixedlayerdepths,-chl_max_depth,
             -hypoxia_depth,-beam_depth,-chl_max,-intergrated_chl,-distance_from_shore,-PAR_1_depth_adj,-day_night_0_1,-density2))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(~PC1, data=sample_dat))
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
    
    # Create the PDF filename string
    save_path="plots/models_by_var/"
    filename <- paste0(save_path,"S_", dat_name, "_PC1_18s", ".png")
    
    # For demonstration purposes, print the filename
    cat("Saving to:", filename, "\n")
    
    
    png(filename, width = 1000, height = 600)
    
    # Your plotting code
    sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
    pp=plot(posterior, par="Lambda", focus.coord = focus, focus.cov = "PC1")+
      labs(title=paste(sizes[counter]),
           x="")
    base::print(pp)
    Sys.sleep(2)
    # Close the PDF device
    dev.off()
    #Save loop plot
    assign(paste0("pp", counter), pp) 
  }
}


all_plot_18s=grid.arrange(pp1,pp2,pp3, nrow=1,
                      bottom = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore")
all_plot_18s



#Now make into a ggplot
#####
#Add size column
posterior_summary_S1 =posterior_summary_S1 %>% mutate(size="0.2-0.5")
posterior_summary_S2 =posterior_summary_S2 %>% mutate(size="0.5-1")
posterior_summary_S3 =posterior_summary_S3 %>% mutate(size="1-2")

posterior_summary_18s=bind_rows(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3)%>%
  filter(covariate != "(Intercept)")

#Plot
all_plot_18s=ggplot(posterior_summary_18s, aes(x = mean, y = coord, color=size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = p2.5, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape=size),size = 8, alpha=0.8) + 
  labs(title="Pibble Model for 18S Taxa vs. PC1",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Size Fraction (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 
all_plot_18s

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q3_taxa_vars/fido_PC1_18s_taxa.pdf"), 
    plot = all_plot_18s,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q3_taxa_vars/fido_PC1_18s_taxa.png"), 
    plot = all_plot_18s,
    width = 16,  # Width in inches
    height = 6  # Height in inches
  )}


#=============#Testing all vars

counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(-"Sample_ID_short")
  meta_vars=names(sample_dat[,1:ncol(sample_dat)])
  
  for (var in meta_vars) {
    formula_string <- paste("~", var)
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
    
    tax <- tax_table(dat)[,c("Genus")] %>% as.data.frame() %>%
      rownames_to_column("Genus2")%>% select(Genus)
    num <- 1:nrow(tax)
    tax <- unname(apply(tax, 1, paste, collapse="_"))
    tax <- paste(tax,sep="_")
    names_categories(posterior) <- tax
    
    
    ##This is the "now what?" part. We have our model, what does it tell us?
    posterior_summary <- summary(posterior, pars="Lambda")$Lambda
    
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
      focus_cov
      
      # Create the PDF filename string
      save_path="plots/models_by_var/"
      filename <- paste0(save_path,"S_", dat_name, "_", var, ".png")
      
      # For demonstration purposes, print the filename
      cat("Saving to:", filename, "\n")
      
      
      png(filename, width = 1000, height = 600)
      
      # Your plotting code
      sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
      pp=plot(posterior, par="Lambda", focus.coord = focus, focus.cov = focus_cov)+
        labs(title=paste(sizes[counter], var))
      base::print(pp)
      Sys.sleep(2)
      # Close the PDF device
      dev.off()
    }
    
    
  }
  
  
}



