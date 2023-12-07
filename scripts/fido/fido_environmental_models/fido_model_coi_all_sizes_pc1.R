library(MicrobeDS)
library(phyloseq)
library(tidyverse)
library(fido)
library(corrplot)
library(here)

here()
#coi read in different sizes
otucoi1=read.csv(here("data/fido/fido_coi_s1_ecdf.csv"), row.names = 1)%>%
  dplyr::select(where(~ !is.na(.[[1]])))%>% 
  select(starts_with("C"))%>%
  rownames_to_column("Hash") %>%
  pivot_longer(-Hash, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Hash) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Hash") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otucoi2=read.csv(here("data/fido/fido_coi_s2_ecdf.csv"), row.names = 1)%>%
  dplyr::select(where(~ !is.na(.[[1]])))%>% 
  select(starts_with("C"))%>%
  rownames_to_column("Hash") %>%
  pivot_longer(-Hash, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Hash) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Hash") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otucoi3=read.csv(here("data/fido/fido_coi_s3_ecdf.csv"), row.names = 1)%>%
  dplyr::select(where(~ !is.na(.[[1]])))%>% 
  select(starts_with("C"))%>%
  rownames_to_column("Hash") %>%
  pivot_longer(-Hash, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Hash) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Hash") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

taxcoi=read.csv(here("data/metazooprunedcoi_tax.csv"),row.names = 1)%>% add_row()
rownames(taxcoi)[nrow(taxcoi)] <- 'other'
#Subset to match 
taxcoi1 = taxcoi %>% filter(rownames(taxcoi) %in% rownames(otucoi1))
taxcoi1=  tax_table(as.matrix(taxcoi1))
taxcoi2 = taxcoi %>% filter(rownames(taxcoi) %in% rownames(otucoi2))
taxcoi2=  tax_table(as.matrix(taxcoi2))
taxcoi3 = taxcoi %>% filter(rownames(taxcoi) %in% rownames(otucoi3))
taxcoi3=  tax_table(as.matrix(taxcoi3))


metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,cycle, max_size)) %>%
  sample_data(.)

dat_1=phyloseq(otucoi1,taxcoi1,metacoi) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otucoi2,taxcoi2,metacoi)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otucoi3,taxcoi3,metacoi)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)


set.seed(899)


counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(-"Sample_ID_short")
  meta_vars=names(sample_dat[,1:ncol(sample_dat)])

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
  
  tax <- tax_table(dat)[,c("Family","Species")]
  hash <- rownames(tax_table(dat))
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,hash,sep="_")
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
    # filter(covariate != "(Intercept)")
  # focus <- posterior_summary[sign(posterior_summary$p25) == sign(posterior_summary$p75),]
  
  focus ##Note that there are zero rows --> no evidence of an effect of potemp2 on any of the taxa
  
  ##This code will only work if there is a sig. result returned
  focus_coord <- unique(focus$coord)
  # focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5)]]
  focus_cov= "PC1"
  focus_cov
  
  # Create the PDF filename string
  save_path="plots/model_multiple_vars/"
  filename <- paste0(save_path,"S_", dat_name, "_multiple_vars_coi.png")
  
  # For demonstration purposes, print the filename
  cat("Saving to:", filename, "\n")
  
  
  png(filename, width = 1000, height = 600)
  
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  pp=plot(posterior, par="Lambda", focus.coord = focus_coord, focus.cov = focus_cov)+
    labs(title=paste(sizes[counter]))
pp
  Sys.sleep(2)
  # Close the PDF device
  dev.off()
  
}


#####
#Add size column
posterior_summary_S1 =posterior_summary_S1 %>% mutate(size="0.2-0.5")
posterior_summary_S2 =posterior_summary_S2 %>% mutate(size="0.5-1")
posterior_summary_S3 =posterior_summary_S3 %>% mutate(size="1-2")

posterior_summary_all=bind_rows(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3)%>%
  filter(covariate != "(Intercept)")

#Plot
ggplot(posterior_summary_all, aes(x = mean, y = coord, color=size)) + 
  geom_point(size = 3) + 
  geom_vline(aes(xintercept = 0), color = "black")+
  geom_errorbarh(aes(xmin = p2.5, xmax = p97.5), height = 0.2) +
  labs(x = "Mean Log-Ratio PC1",y="", color = "Size Class") +
  theme_minimal()+
  facet_wrap(~size)


