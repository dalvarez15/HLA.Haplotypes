############################################################
## Estimation of HLA-DRB1~DQB1 haplotypes in individual   ##
## samples using NMDP haplotype frequencies               ##
############################################################

## Call estimate_haplotypes() function for HLA-DRB1~DQB1 in the whole panel using apply()
DQDR.estimated.haplotypes <- data.frame()
DQDR.estimated.haplotypes <- apply(df.2fields.typings[which(df.2fields.typings$`HLA-DRB1`!= "NA" & df.2fields.typings$`HLA-DQB1` != "NA"),], 
         1, estimate_haplotypes, "HLA-DRB1", "HLA-DQB1", nmdp.DQB1DRB1, DQDR.estimated.haplotypes)
DQDR.estimated.haplotypes <- do.call(rbind, DQDR.estimated.haplotypes)

# Clean output
colnames(DQDR.estimated.haplotypes)[c(1,2)] <- c("haplotype.y", "haplotype.x")
DQDR.estimated.haplotypes <- DQDR.estimated.haplotypes %>%
  select(sample, `mother/child`, haplotype.y, haplotype.x, gene1.alleles.x, gene2.alleles.x, CAU_freq.x, CAU_rank.x,
         LD.x, "D'.x", allele1prop.x, allele2prop.x, gene1.alleles.y, gene2.alleles.y, CAU_freq.y, CAU_rank.y, 
         LD.y, "D'.y", allele1prop.y, allele2prop.y, metric, normalized, unique.alleles)

# Select haplotype combinations with best metrics per sample
DQDR.estimated.haplotypes <- DQDR.estimated.haplotypes[which(DQDR.estimated.haplotypes$normalized == 1),]
# HLA-DQB1~DRB1 haplotypes estimation output
head(BC.estimated.haplotypes)

## Allele exchanging method in search of more frequent haplotypes
# Call allele_exchanging_haplotypes() function in all estimated haplotypes using
# apply()
DQDR.exchanged.haplotypes <- data.frame() 
DQDR.exchanged.haplotypes <- apply(DQDR.estimated.haplotypes, 1, allele_exchanging_haplotypes, nmdp.DQB1DRB1, DQDR.exchanged.haplotypes)
DQDR.exchanged.haplotypes <- do.call(rbind, DQDR.exchanged.haplotypes)

# Clean output
DQDR.exchanged.haplotypes <- DQDR.exchanged.haplotypes %>%
  select("sample", "mother/child", "DRB1", "DQB1", "CAU_freq", "CAU_rank", "prev.haplotype",
         "prev.allele", "prev.CAU_rank", "unique.alleles")
# Allele exchanging method output for complete HLA-DQB1~DRB1 haplotypes
head(DQDR.exchanged.haplotypes)

## Allele exchanging method to exchanging only subtypes (second field of alleles)
# Call allele_exchanging_subtypes() function in all estimtated haplotypes using
# apply()
DQDR.exchanged.subtypes <- data.frame(matrix(ncol=24, nrow=0))
DQDR.exchanged.subtypes <- apply(DQDR.estimated.haplotypes, 1, allele_exchanging_subtypes, nmdp.DQB1DRB1, DQDR.exchanged.subtypes)
DQDR.exchanged.subtypes <- do.call(rbind, DQDR.exchanged.subtypes)
# Fix output
DQDR.exchanged.subtypes <- DQDR.exchanged.subtypes %>%
  select("sample", "mother/child","DQB1", "DRB1", "CAU_freq", "CAU_rank", "LD", "D'", allele1prop, allele2prop, prev.haplotype,
         "prev.allele", "prev.CAU_freq", "prev.CAU_rank", prev.LD, "prev.D'", "unique.alleles")
# Allele exchanging method output for HLA-DQB1~DRB1 subtypes
head(DQDR.exchanged.subtypes)

# Create individual sample.id to obtain 
# total number of samples in which more frequent haplotypes were found
DQDR.estimated.haplotypes$sample.id <- paste0(DQDR.estimated.haplotypes$sample, DQDR.estimated.haplotypes$`mother/child`)
DQDR.exchanged.haplotypes$sample.id <- paste0(DQDR.exchanged.haplotypes$sample, DQDR.exchanged.haplotypes$`mother/child`)
DQDR.exchanged.subtypes$sample.id <- paste0(DQDR.exchanged.subtypes$sample, DQDR.exchanged.subtypes$`mother/child`)

# Plot results
pie.DQDR.exchanged.haplotypes <- data.frame(
  group =  c("Positive", "Negative"),
  value = c(length(unique(DQDR.exchanged.haplotypes$sample.id)), 
            length(unique(DQDR.estimated.haplotypes$sample.id))-length(unique(DQDR.exchanged.haplotypes$sample.id)))
)
pie.DQDR.exchanged.haplotypes <- ggplot(pie.DQDR.exchanged.haplotypes, aes(x="", y=value, fill=reorder(group, -value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(fill = "Reported more frequent haplotypes",
       x = NULL,
       y = NULL,
       title = "HLA-DQB1~DRB1 alleles")+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))
pie.DQDR.exchanged.haplotypes

pie.DQDR.exchanged.subtypes <- data.frame(
  group =  c("Positive", "Negative"),
  value = c(length(unique(DQDR.exchanged.subtypes$sample.id)), 
            length(unique(DQDR.estimated.haplotypes$sample.id))-length(unique(DQDR.exchanged.subtypes$sample.id)))
)
pie.DQDR.exchanged.subtypes <- ggplot(pie.DQDR.exchanged.subtypes, aes(x="", y=value, fill=reorder(group, value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(fill = "Reported more frequent haplotypes",
       x = NULL,
       y = NULL,
       title = "HLA-DQB1~DRB1 subtypes")+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))
pie.DQDR.exchanged.subtypes

# Merge allele exchanging plot results from HLA-B~C, HLA-DQB1~DRB1 (complete alleles and subtypes) in one plot
grid_arrange_shared_legend(pie.DQDR.exchanged.haplotypes, pie.DQDR.exchanged.subtypes, pie.BC.exchanged.haplotypes, pie.BC.exchanged.subtypes,
                           nrow = 2,  
                           ncol = 2,
                           top = textGrob("Allele exchanging results per sample", gp=gpar(fontsize=16), vjust = 0.5))
