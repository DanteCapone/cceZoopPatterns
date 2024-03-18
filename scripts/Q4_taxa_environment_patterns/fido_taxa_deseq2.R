
#Deseq2 comparison for taxonomic subsets identified by cluster comparison

#For 18s, the 2 clusters lined up best for all, size 1 and 2

#Load packages
librarian::shelf(tidyverse, googledrive, stringr,here,DESeq2,phyloseq,
                 extrafont, vegan)
dev.off()
#Load the data
#18s reads
otu18s=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_otu.csv")) %>%
  column_to_rownames("Hash")%>%
  dplyr::select(where(~ !is.na(.[[1]]))) %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

#Add pseudocount
pseudocount <- 1
otu18s[otu18s == 0] <- pseudocount


tax18s=read.csv(here("data/phyloseq_bio_data/18s/metazoopruned18s_tax.csv"))
tax18s=tax_table(as.matrix(tax18s %>% column_to_rownames("Hash")))
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(clust_group,cycle)) %>%
  sample_data(.)

#

zhan=phyloseq(otu18s,tax18s,meta18s)
zhan_05=zhan %>% subset_samples(.,max_size==0.5)
zhan_1=zhan %>% subset_samples(.,max_size==1)

#deseq2
phy_sel=zhan_1
dds = phyloseq_to_deseq2(phy_sel, ~offshore_onshore)
dds = DESeq(dds, test="Wald", fitType="parametric")





#########PLot

#Normalize

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
counts(dds, normalized=TRUE)

# Extract differential abundance results
res <- results(dds)



#
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy_sel)[rownames(sigtab), ], "matrix"))
head(sigtab)

#Convert names
sigtab =sigtab %>%
  mutate(Species = if_else(Species == "Species", "Unidentified Species", Species))


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))


#filter
sigtab_p=sigtab %>% filter(Family != "NA")%>% filter(Species != "NA") %>% filter(Order != "NA")

min(levels(sigtab_p$Family))

#Order
ggplot(sigtab_p, aes(x=Family, y=log2FoldChange, color=Order))+
  geom_hline(yintercept = 0, color = "black", size = 1, alpha=0.7) + geom_point(size=6)+ coord_flip()+
  theme_classic()+
  labs(title=c(paste(sample_data(phy_sel)$Sizefractionmm[1],"mm zhan")))+
geom_text(aes(y = -1, x = min(levels(sigtab_p$Family)), label = "Offshore"), hjust = 1, vjust = 0.5, color = "black", fontface = "bold") +
  geom_text(aes(y = 1,  x= min(levels(sigtab_p$Family)), label = "Onshore"), hjust = 0, vjust = 0.5, color = "black", fontface = "bold")




