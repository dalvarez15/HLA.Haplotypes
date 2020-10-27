############################################################
## Estimation of HLA-B~C haplotypes in individual samples ##
## using NMDP haplotype frequencies                       ##
############################################################

## Call estimate_haplotypes() function for HLA-B~C in the whole panel using apply()
BC.estimated.haplotypes <- data.frame()
BC.estimated.haplotypes <- apply(df.2fields.typings[which(df.2fields.typings$`HLA-C`!= "NA" & df.2fields.typings$`HLA-C` != "NA"),], 
          1, estimate_haplotypes, "HLA-C", "HLA-B", nmdp.BC, BC.estimated.haplotypes)
BC.estimated.haplotypes <- do.call(rbind, BC.estimated.haplotypes)

# Clean output
colnames(BC.estimated.haplotypes)[c(1,2)] <- c("haplotype.y", "haplotype.x")
BC.estimated.haplotypes <- BC.estimated.haplotypes %>%
  select(sample, `mother/child`, haplotype.y, haplotype.x, gene1.alleles.x, gene2.alleles.x, CAU_freq.x, CAU_rank.x,
         LD.x, "D'.x", allele1prop.x, allele2prop.x, gene1.alleles.y, gene2.alleles.y, CAU_freq.y, CAU_rank.y, 
         LD.y, "D'.y", allele1prop.y, allele2prop.y, metric, normalized, unique.alleles)

# Select haplotype combinations with best metrics per sample
BC.estimated.haplotypes <- BC.estimated.haplotypes[which(BC.estimated.haplotypes$normalized == 1),]
# HLA-B~C haplotypes estimation output
head(BC.estimated.haplotypes)

## Allele exchanging method in search of more frequent haplotypes
# Call allele_exchanging_haplotypes() function in all estimated haplotypes using
# apply()
BC.exchanged.haplotypes <- data.frame() 
BC.exchanged.haplotypes <- apply(BC.estimated.haplotypes, 1, allele_exchanging_haplotypes, nmdp.BC, BC.exchanged.haplotypes)
BC.exchanged.haplotypes <- do.call(rbind, BC.exchanged.haplotypes)
  
# Clean output
BC.exchanged.haplotypes <- BC.exchanged.haplotypes %>%
  select("sample", "mother/child", "C", "B", "CAU_freq", "CAU_rank", "prev.haplotype",
         "prev.allele", "prev.CAU_rank", "unique.alleles")
# Allele exchanging method output for complete HLA-B~C haplotypes
head(BC.exchanged.haplotypes)

## Allele exchanging method to exchanging only subtypes (second field of alleles)
# Call allele_exchanging_subtypes() function in all estimated haplotypes using
# apply()
BC.exchanged.subtypes <- data.frame()
BC.exchanged.subtypes <- apply(BC.estimated.haplotypes,1,allele_exchanging_subtypes, nmdp.BC, BC.exchanged.subtypes)
BC.exchanged.subtypes <- do.call(rbind, BC.exchanged.subtypes)
# Clean output
BC.exchanged.subtypes <- BC.exchanged.subtypes %>%
  select("sample", "mother/child","C", "B", "CAU_freq", "CAU_rank", "LD", "D'", allele1prop, allele2prop, prev.haplotype,
         "prev.allele", "prev.CAU_freq", "prev.CAU_rank", prev.LD, "prev.D'","unique.alleles")
# Allele exchanging method output for HLA-B~C subtypes
head(BC.exchanged.subtypes)

# Create individual sample.id to obtain 
# total number of samples in which more frequent haplotypes were found
BC.estimated.haplotypes$sample.id <- paste0(BC.estimated.haplotypes$sample, BC.estimated.haplotypes$`mother/child`)
BC.exchanged.haplotypes$sample.id <- paste0(BC.exchanged.haplotypes$sample, BC.exchanged.haplotypes$`mother/child`)
BC.exchanged.subtypes$sample.id <- paste0(BC.exchanged.subtypes$sample, BC.exchanged.subtypes$`mother/child`)

# Plot results
pie.BC.exchanged.haplotypes <- data.frame(
  group = c("Positive", "Negative"),
  value = c(length(unique(BC.exchanged.haplotypes$sample.id)), 
                   length(unique(BC.estimated.haplotypes$sample.id))-length(unique(BC.exchanged.haplotypes$sample.id)))
)
pie.BC.exchanged.haplotypes <- ggplot(pie.BC.exchanged.haplotypes, aes(x="", y=value, fill=reorder(group, -value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(fill = "",
       x = NULL,
       y = NULL,
       title = "HLA-B~C alleles")+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))
pie.BC.exchanged.haplotypes

pie.BC.exchanged.subtypes <- data.frame(
  group = c("Positive", "Negative"),
  value = c(length(unique(BC.exchanged.subtypes$sample.id)), 
                   length(unique(BC.estimated.haplotypes$sample.id))-length(unique(BC.exchanged.subtypes$sample.id)))
)
pie.BC.exchanged.subtypes <- ggplot(pie.BC.exchanged.subtypes, aes(x="", y=value, fill=reorder(group, value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(fill = "Reported more frequent haplotypes",
       x = NULL,
       y = NULL,
       title = "HLA-B~C subtypes")+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))
pie.BC.exchanged.subtypes

# Go to script "4_DQB1DRB1_estimation.R" to find plot with all allele exchanging results
# for both analysed haplotypes, search of complete alleles and restriction to allele subtypes
