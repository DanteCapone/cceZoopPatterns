#Compare amp effs with abundance in size


amp_effs_18s=read.csv(here("data/amp_effs/all_amp_effs_18s.csv"))%>%
  mutate(Family = str_remove(Lambda.coord, "clr_"))



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
  mutate(Sample_ID=Sample_ID_dot)


fido_18s_raw=fido_18s_merged_raw %>% 
  rownames_to_column("Family")%>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>% 
  left_join(env_metadata, by="Sample_ID") %>%
  group_by(Sample_ID) %>% 
  mutate(total_reads=sum(n_reads)) %>% 
  ungroup() %>% 
  mutate(reads_proportions=n_reads/total_reads) %>% 
  group_by(Family, size_fraction) %>% 
  summarise(n_reads=sum(reads_proportions)) 


amp_effs_and_raw=amp_effs_18s %>% 
  merge(fido_18s_raw, by=c("Family","size_fraction"))


amp_effs_and_raw %>% 
  filter(n_reads<4) %>% 
  ggplot(aes(x=n_reads, y=Lambda.mean,color=Family))+
  geom_point(aes())+
  facet_wrap(~size_fraction, nrow=3)+
  stat_cor(method="pearson", label.x = 5, label.y = -0.15,size=3)+
  theme_classic()

