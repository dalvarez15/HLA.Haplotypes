##############################################
# Analysis of allele dropout data and        #
# samples with at least one real homozygous  #
##############################################

# Load small panel with allele dropout data
allele.dropout.df <- read_excel("Input/AlleleDropOutData.xlsx")

# Curate data.frame
allele.dropout.df$`HLA-B` <- gsub('G', '', allele.dropout.df$`HLA-B`)
allele.dropout.df$`HLA-B 2` <- gsub('G', '', allele.dropout.df$`HLA-B 2`)
allele.dropout.df$`HLA-C` <- gsub('G', '', allele.dropout.df$`HLA-C`)
allele.dropout.df$`HLA-C 2` <- gsub('G', '', allele.dropout.df$`HLA-C 2`)

# Reduce typing resolution to two fields
allele.dropout.df[,c(2,3,4,5)]  <- reduce_typing_resolution(allele.dropout.df[,c(2,3,4,5)], 2)

# Prepare data.frame for translation to two fields plus g group
colnames(allele.dropout.df) <- str_replace(colnames(allele.dropout.df), "HLA-", "")
colnames(allele.dropout.df) <- str_replace(colnames(allele.dropout.df), " 2", "")
colnames(allele.dropout.df)[seq(2,5)] <- paste0(colnames(allele.dropout.df)[seq(2,5)], "*")
allele.dropout.df[,seq(2,5)] <- Map(paste0, names(allele.dropout.df[,seq(2,5)]), allele.dropout.df[,seq(2,5)])
# Translate into g group calling translate_g_group() with apply()
allele.dropout.df[,seq(2,5)] <- apply(allele.dropout.df[,seq(2,5)], c(1,2), translate_g_group, translation) 
colnames(allele.dropout.df)[seq(1,5)] <- c("no.","HLA-B", "HLA-B 2", "HLA-C", "HLA-C 2")

# Estimate haplotypes in the data.frame
AlDO.haplotypes <- data.frame(matrix(ncol=37, nrow=0))
AlDO.haplotypes <- apply(allele.dropout.df, 1, estimate_haplotypes, "HLA-C", "HLA-B", nmdp.BC, AlDO.haplotypes)
AlDO.haplotypes <- do.call(rbind, AlDO.haplotypes)

# Clean output
colnames(AlDO.haplotypes)[c(1,2)] <- c("haplotype.y", "haplotype.x")
AlDO.haplotypes <- AlDO.haplotypes %>%
  select(haplotype.y, haplotype.x, gene1.alleles.x, gene2.alleles.x, CAU_freq.x, CAU_rank.x,
         LD.x, "D'.x", gene1.alleles.y, gene2.alleles.y, CAU_freq.y, CAU_rank.y, 
         LD.y, "D'.y",  metric, normalized, sample, unique.alleles)

# Call allele_exchanging_homozygous() using apply()
AlDO.experiment <- data.frame()
AlDO.experiment <- apply(AlDO.haplotypes, 1, allele_exchanging_homozygous, nmdp.BC, AlDO.experiment)
AlDO.experiment <- do.call(rbind, AlDO.experiment)

# Output of allele exchanging method for samples with at least one homzoygous 
head(AlDO.experiment)

# Add unique sampleIDs
AlDO.experiment$sampleID <- paste0(AlDO.experiment$sample, AlDO.experiment$prev.haplotype,
                                   AlDO.experiment$prev.allele, AlDO.experiment$propsum)
# Convert into data.table
AlDO.experiment <- data.table(AlDO.experiment)
# Obtain propsum metric
propsums.DO <- unique(AlDO.experiment, by="sampleID")
propsums.DO$gene <- as.factor("HLA-B~C Allele Dropouts")

## Analyisis of real homozygous data with allele_exchanging_homozygous() method
# In HLA-B~C
homo.BC.experiment <- data.frame()
homo.BC.experiment <- apply(BC.estimated.haplotypes[which(BC.estimated.haplotypes$unique.alleles<4),],1,allele_exchanging_homozygous, nmdp.BC, homo.BC.experiment)
homo.BC.experiment <- do.call(rbind, homo.BC.experiment)
homo.BC.experiment$sampleID <- paste0(homo.BC.experiment$sample, homo.BC.experiment$`mother/child`, homo.BC.experiment$prev.haplotype,
                                      homo.BC.experiment$prev.allele, homo.BC.experiment$propsum)
homo.BC.experiment <- data.table(homo.BC.experiment)
propsums.homo.BC <- unique(homo.BC.experiment, by="sampleID")
propsums.homo.BC$gene <- as.factor("HLA-B~C Homozygous")

# In HLA-DQB1~DRB1
homo.DQDR.experiment <- data.frame()
homo.DQDR.experiment <- apply(DQDR.estimated.haplotypes[which(DQDR.estimated.haplotypes$unique.alleles<4),],1,allele_exchanging_homozygous, nmdp.DQB1DRB1, homo.DQDR.experiment)
homo.DQDR.experiment <- do.call(rbind, homo.DQDR.experiment)
homo.DQDR.experiment$sampleID <- paste0(homo.DQDR.experiment$sample, homo.DQDR.experiment$`mother/child`, 
                                        homo.DQDR.experiment$prev.haplotype,
                                        homo.DQDR.experiment$prev.allele, homo.DQDR.experiment$propsum)
homo.DQDR.experiment <- data.table(homo.DQDR.experiment)
propsums.homo.DQDR <- unique(homo.DQDR.experiment, by="sampleID")
propsums.homo.DQDR$gene <- as.factor("HLA-DQB1~DRB1 Homozygous")

## Create boxplot to compare obtained propsums metric in allele dropout data, 
#  B-C and DQDR real homozygous
dropout.vs.homo <- rbind(propsums.DO[,c("propsum", "gene")], 
                         propsums.homo.BC[,c("propsum", "gene")],
                         propsums.homo.DQDR[,c("propsum", "gene")])

boxplot.dropouts <- ggplot(dropout.vs.homo, aes(x=gene, y=propsum)) + 
  geom_boxplot()+
  labs(x="Haplotype and category",
       y="Propsum (%)",
       title="Propsum metrics from homozygous allele exchanging method"
  )
boxplot.dropouts
