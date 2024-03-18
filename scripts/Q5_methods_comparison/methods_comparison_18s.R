


#Plot PCR-Bias mitigated data as a function of PC1
librarian::shelf(tidyverse, googledrive, stringr,here,gridextra,phyloseq,
                 extrafont, RColorBrewer)

# Packages and Functions --------------------------------------------------




#Add functions for myself
source(("scripts/helpful_functions/treemap_funs_Capone.R"))
source("scripts/helpful_functions/phyloseq_mapping_funs.R")
source("scripts/helpful_functions/general_helper_functions.R")



# Loading in the data -----------------------------------------------------


saving=0
#Metadata
##Lat and lon are funky so add back in from metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  select(-X, -Sizefractionmm,max_size) %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)

#Add depth
depths=read.csv(here("data/physical_environmental_data/sample_depths.csv")) %>%
  select(-X) %>%
  mutate(Sample_ID_short=Sample_ID)


#Volume filtered (add to metadata)
volume_filtered=read.csv(here("data/raw_data/biomass/p2107_bt_volume_filtered.csv"))

#Dryweights
dryweights=read.csv("data/raw_data/biomass/dryweights_forzoopmetab.csv") %>%
  mutate(biomass_dry=8/3*biomass_dry)

env_metadata=metadata %>% left_join(.,volume_filtered, by="Sample_ID_short") %>%
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x)





# Methods Comparison ------------------------------------------------------

#Taxa
zhan_taxa=read.csv(here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv")) %>%
  column_to_rownames("Family") %>% 
  mutate(Hash=X) %>%
  select(-X, -Species, -Genus)

#Join PCR and RRA dataframes

#===== PCR Bias Mitigated Proportion Data
#Predicted proportions
fido_s1=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s1_phy.csv")) %>%
  select(-X)
fido_s2=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s2_phy.csv")) %>%
  select(-X)
fido_s3=read.csv(here("data/predicted_og/predicted_og_18s_02_22_2024_s3_phy.csv")) %>%
  select(-X)
#Merge
final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 

#Make final dataframe
phy_taxa_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr=phy_taxa_pcr %>% mutate(Family=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Family"), by="Family") %>%
  filter(Order=="Calanoida")

pcr_join=taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction,PC1,cycle) %>%
  summarise(n_reads_pcr=sum(n_reads))


# =========== Raw Reads Realtive Abundance Data using taxa that went into fido model




#Predicted proportions
fido_s1_raw=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv")) %>% 
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
env_metadata_phy=env_metadata %>%
  column_to_rownames("Sample_ID_dot")

#Make phyloseq objects

#18s
# OTU = otu_table(as.matrix(zhan_otu), taxa_are_rows = TRUE)
OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(env_metadata_phy)

#USe proportions
phy_18s=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>%
  select(-asv_code)

#filter to calanoid copepods
taxa_sel="Calanoida"

taxa_raw=phy_18s %>% filter(Order==taxa_sel)

taxa_raw %>%
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)->pcr_and_raw

#Add difference column
pcr_and_raw_18s=pcr_and_raw %>%
  mutate(difference_pcr_raw=n_reads_pcr-n_reads_raw)

#Plot relative abundance differences
pcr_and_raw_18s %>%
  # filter(Species %in% phy_spp$taxa) %>%
  ggplot(aes(x = as.factor(PC1), y = difference_pcr_raw, fill = cycle)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  labs(title = "PCR-Raw",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance (Corrected-Raw)",
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  # coord_cartesian(ylim=c(-0.5,0.8))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)-> pcr_raw_plot
pcr_raw_plot

if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/pcr_rra_relative_abundance_diff.pdf"), 
    plot = pcr_raw_plot,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )}


# Zooscan -----------------------------------------------------------------


#===== Zooscan comparison

#Need to modify string category for joining
size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)


#Read in processed Zooscan data and look at biomass
zooscan_by_sample=read.csv(here("data/Zooscan/zooscan_by_sample_biomass.csv")) %>%
  select(-X) %>%
  mutate(Sample_ID=sample_id) %>%
  mutate(size_fraction = case_when(
    size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
    TRUE ~ NA_real_)) %>%
  group_by(Sample_ID)

#Add biomass sum, calanoid biomass and proportion of calanoid biomass
taxa_sel="Calanoida"
zooscan_taxa=zooscan_by_sample %>%
  filter(object_annotation_category=="Calanoida") 


## ==== Biomass & Biomass proportions plots === #
labels_for_map=biomass_map %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

zooscan_taxa %>%
  filter(size_fraction!=5) %>%
  ggplot(aes(x = as.factor(PC1), y = biomass_prop_taxa, fill = cycle)) +
  geom_bar(stat = "identity", width=0.8) +
  labs(title = "Calanoid Copepod Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->taxa_biomass_prop_zooscan
taxa_biomass_prop_zooscan

#Save
if (saving==1) {
# ggsave(
#   filename = here("plots/methods_comparison/zooscan_calanoid_biomass_proportions.pdf"), 
#   plot = taxa_biomass_zooscan,
#   width = 8,  # Width in inches
#   height = 6  # Height in inches
# )
}


#Biomass
custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  
zooscan_taxa %>%
  filter(size_fraction!=5) %>%
  ggplot(aes(x = as.factor(PC1), y = dryweight_C_mg_m2_taxa, fill = cycle)) +
  geom_bar(stat = "identity", width=0.8) +
  labs(title = "Calanoid Copepod Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->taxa_biomass_zooscan
taxa_biomass_zooscan

#Save
if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/zooscan_calanoid_biomass_proportions.pdf"), 
  plot = taxa_biomass_zooscan,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)
}

#========== COMPARE: Make combined dataframe for comparing all 3 methods
pcr_raw_zoo_18s=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_18s, by=c("PC1","size_fraction"))%>%
  #Compare absolute biomasses from Drymass
  mutate(difference_pcr_zoo_biomass=n_reads_pcr*biomass_mg_m2*0.37-dryweight_C_mg_m2_taxa,
         difference_raw_zoo_biomass=n_reads_raw*biomass_mg_m2*0.37-dryweight_C_mg_m2_taxa) %>%
  #Compare proportions using reads
  mutate(difference_pcr_zoo_biomass_prop=n_reads_pcr-biomass_prop_taxa,
         difference_raw_zoo_biomass_prop=n_reads_raw-biomass_prop_taxa)







#=============Grouped Barplot site abundances
#Format Long
pcr_raw_zoo_18s_long <- pivot_longer(pcr_raw_zoo_18s, 
                                     cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05*relative_abundance, "xx",NA))




# ANOVA test
anova_result <- aov(relative_abundance ~ Method * as.factor(PC1), data = pcr_raw_zoo_18s_long)
anova_summary <- summary(anova_result)
print(anova_summary)

# Kruskal-Wallis test (non-parametric alternative)
kruskal_result <- pcr_raw_zoo_18s_long %>%
  group_by(Sample_ID.y) %>%
  do(kruskal_test = kruskal.test(relative_abundance ~ Method, data = .))

# Print results
print(kruskal_result$kruskal_test)

# Perform Tukey's HSD post hoc test
tukey_result <- TukeyHSD(anova_result)

# Print the results
print(tukey_result)

#Plot grouped bar plot
pcr_raw_zoo_18s_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_point(position = position_dodge(width = 0.9), aes(shape = Method), size = 3) +  # Add points
  # geom_line(position = position_dodge(width = 0.9), aes(group = Method, color = Method), size = 2) +  # Add lines with colors
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Biomass ",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#70BF41", "#4F86F7", "#F78D4F"),
                    labels = c("Zooscan", "PCR-corrected", "Raw Reads")) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+# Define linetypes for each Method
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F")) -> grouped_bar_all  # Define colors for each Method

grouped_bar_all

if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s.pdf"),
  # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
  plot = grouped_bar_all,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s.png"),
  # filename = here("plots/methods_comparison/line_relative_abundance_diff.png"),
  plot = grouped_bar_all,
  width = 8,  # Width in inches
  height = 6  # Height in inches
)
}

## Correlation plot

# Check distribution
pcr_raw_zoo_18s %>%
  ggplot(., aes(x = asin(sqrt(n_reads_pcr)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

pcr_raw_zoo_18s %>%
  ggplot(., aes(x = asin(sqrt(n_reads_raw)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

pcr_raw_zoo_18s %>%
  ggplot(., aes(x = asin(sqrt(biomass_prop_taxa)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels



# Correlation


#ADd clusters to df for plotting groups
clusters=read.csv(here("data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_18s=pcr_raw_zoo_18s %>%
  left_join(.,clusters,by="PC1")
#Check outliers
# Calculate quantiles
quantiles <- quantile(pcr_raw_zoo_18s$biomass_prop_taxa, c(0.25, 0.75))
IQR <- quantiles[2] - quantiles[1]

# Define lower and upper bounds (e.g., using 1.5*IQR)
lower_bound <- quantiles[1] - 1.5 * IQR
upper_bound <- quantiles[2] + 1.5 * IQR
pcr_raw_zoo_18s

pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "PCR Bias-Mitigated Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.3)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr
zoo_vs_pcr

if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_size.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_cycle.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clust.pdf"),
  # filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_sig_diffs.pdf"), 
  plot = zoo_vs_pcr,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_size.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_cycle.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clust.png"),
  plot = zoo_vs_pcr,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)
}

#### =====  
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_size.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_cycle.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clust.pdf"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_cycle.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_size.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_cycle.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clust.png"),
  # filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clust.png"),
  
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)
}


### PCR vs Raw
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=n_reads_pcr, y=n_reads_raw))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.2)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw
ggsave(
  filename = here("plots/methods_comparison/pcr_vs_raw_correlation.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/pcr_vs_raw_correlation.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)



# ANCOVA ------------------------------------------------------------------



## ANCOVA with groupings
library(rstatix)
#18s raw-RA
ancova_result_18s_raw <- pcr_raw_zoo_18s %>%
  filter(!is.na(n_reads_raw))%>%
  mutate(t1 = ifelse(cycle.y == "T1", cycle.y, "other")) %>%
  lm(asin(sqrt(n_reads_raw)) ~ asin(sqrt(biomass_prop_taxa))+size_fraction, data = .)

#Can remove offshore_onshore and interactions for cycle
anova(ancova_result_18s_raw)
# Diagnostic plots
par(mfrow=c(2,2)) # Create a 2x2 layout for the plots
plot(ancova_result_18s_raw) # Plot diagnostic plots

#18s PCR-RA
ancova_result_18s_pcr <- pcr_raw_zoo_18s %>%
  filter(!is.na(n_reads_raw))%>%
  mutate(t1 = ifelse(cycle.y == "T1", cycle.y, "other")) %>%
  lm(asin(sqrt(n_reads_pcr)) ~ asin(sqrt(biomass_prop_taxa))*size_fraction*cycle.y*offshore_onshore, data = .)

#Can remove offshore_onshore and interactions for cycle
anova(ancova_result_18s_pcr)
# Diagnostic plots
par(mfrow=c(2,2)) # Create a 2x2 layout for the plots
plot(ancova_result_18s_pcr) # Plot diagnostic plots





