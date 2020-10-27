#########################################################################################
# Infer haplotypes in the whole panel at 4 fields and 2 fields + 'g' group resolutions  #
#########################################################################################

## Get vectors with HLA genes names for variables' names 
genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB345", 
           "HLA-DQB1", "HLA-DQA1", "HLA-DPB1", "HLA-DPA1", "HLA-G")
genes2 <- c("HLA-A", "HLA-A 2", "HLA-B", "HLA-B 2", "HLA-C", "HLA-C 2", 
            "HLA-DRB1", "HLA-DRB1 2", "HLA-DRB345", "HLA-DRB345 2", "HLA-DQB1", 
            "HLA-DQB1 2", "HLA-DQA1", "HLA-DQA1 2", "HLA-DPB1", "HLA-DPB1 2", 
            "HLA-DPA1", "HLA-DPA1 2", "HLA-G", "HLA-G 2")

###########################################
## 4 fields / high resolution haplotypes ##
###########################################

# Prepare input for infer_haplotypes() function 
df.4fields.typings <- df.sample.overview[genes2]
colnames(df.4fields.typings) <- genes2
colnames(df.4fields.typings) <- str_replace(colnames(df.4fields.typings), "HLA-", "")
colnames(df.4fields.typings) <- str_replace(colnames(df.4fields.typings), " 2", "")
colnames(df.4fields.typings) <- str_replace(colnames(df.4fields.typings), "345", "")
colnames(df.4fields.typings)[-c(9,10)] <- paste0(colnames(df.4fields.typings)[-c(9,10)], "*")
df.4fields.typings[] <- Map(paste0, names(df.4fields.typings), df.4fields.typings)
colnames(df.4fields.typings) <- genes2
df.4fields.typings <- cbind(df.sample.overview[,c(2,3)], df.4fields.typings)
df.4fields.typings[] <- lapply(df.4fields.typings, function(x) replace(x, grepl("NA", x), "NA"))

# Separate mothers and children typings in two data.frames 
df.4fields.typings.mothers <- df.4fields.typings[which(df.4fields.typings$`mother/child`=="m"),]
df.4fields.typings.children <- df.4fields.typings[which(df.4fields.typings$`mother/child`=="c"),]
# Call infer_haplotypes() using apply() for all samples
df.4fields.haplotypes <- data.frame(matrix(ncol=13, nrow=0))
colnames(df.4fields.haplotypes) <- c("Haplotype", "Sample No.", genes, "Haplotype Length")
df.4fields.haplotypes <- apply(df.4fields.typings.children, 1, infer_haplotypes, df.4fields.typings.mothers,
                               df.4fields.haplotypes)
df.4fields.haplotypes <- do.call(rbind, df.4fields.haplotypes)
#Full resolution (4 fields) haplotypes in the whole panel
head(df.4fields.haplotypes)

# Calculate haplotype frequencies and LD metrics in 4 fields inferred haplotypes 
BC.4fields.frequencies <- get_frequencies(df.4fields.haplotypes, gene1 = "HLA-B", gene2="HLA-C")
DQB1DRB1.4fields.frequencies <- get_frequencies(df.4fields.haplotypes, gene1 = "HLA-DRB1", gene2="HLA-DQB1")

# Get pie chart showing number of categories assigned to all loci (9120 in total)
non_unambiguous_loci <- (9120 - sum(df.4fields.haplotypes[,genes]=="Discordant",
                                     df.4fields.haplotypes[,genes]=="Ambiguous",
                                     df.4fields.haplotypes[,genes]=="Unknown",
                                     df.4fields.haplotypes[,genes]=="NA"))
pie.inferred.haplotypes <- data.frame(
  group = c("Allele unambiguously assigned", "Unknown & NA", "Ambiguous", "Discordant"),
  value = c(non_unambiguous_loci, sum(df.4fields.haplotypes[,genes]=="Unknown", 
            df.4fields.haplotypes[,genes]=="NA"), sum(df.4fields.haplotypes[,genes]=="Ambiguous"), 
            sum(df.4fields.haplotypes[,genes]=="Discordant"))
)
pie.inferred.haplotypes$percentage <- round(100 * pie.inferred.haplotypes$value 
                                            / sum(pie.inferred.haplotypes$value), digits = 1)
pie.inferred.haplotypes<- ggplot(pie.inferred.haplotypes, aes(x="", y=percentage, fill=reorder(group, -percentage)))+
  geom_bar(width = 0.5, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  geom_text(aes(label = scales::percent(value / sum(value)),),
            position = position_stack(vjust = 0.5))+
  labs(fill = NULL,
       x = NULL,
       y = NULL,
       title = "Haplotype categories assigned to all loci")+
  theme_minimal()
pie.inferred.haplotypes # show pie chart

#########################################################
## 2 fields + 'g' group haplotypes (NMDP nomenclature) ##
#########################################################

# Reduce typing resolution to 2 fields
df.2fields.typings <- reduce_typing_resolution(df.sample.overview[genes2], 2)
df.2fields.typings <- as.data.frame(matrix(unlist(df.2fields.typings), 
                                   nrow=length(unlist(df.2fields.typings[1]))))
colnames(df.2fields.typings) <- genes2
colnames(df.2fields.typings) <- str_replace(colnames(df.2fields.typings), "HLA-", "")
colnames(df.2fields.typings) <- str_replace(colnames(df.2fields.typings), " 2", "")
colnames(df.2fields.typings) <- str_replace(colnames(df.2fields.typings), "345", "")
colnames(df.2fields.typings)[-c(9,10)] <- paste0(colnames(df.2fields.typings)[-c(9,10)], "*")

# Translate into two fields + g group nomenclature 
# Import translation to two fields G groups table 
translation <- read_excel("Input/NMDP_translation.xls")
translation[translation=="DRB1*13:01g"] <- "DRB1*13:01"
translation <- data.table(translation)
# Translate into g group
df.2fields.typings[] <- Map(paste0, names(df.2fields.typings), df.2fields.typings)
df.2fields.typings <- apply(df.2fields.typings, c(1,2), translate_g_group, translation) 
df.2fields.typings <- as.data.frame(df.2fields.typings)
colnames(df.2fields.typings) <- genes2
df.2fields.typings <- cbind(df.sample.overview[,c(2,3)], df.2fields.typings)

# Separate mothers and children typings into two data.frames 
df.2fields.mothers <- df.2fields.typings[which(df.2fields.typings$`mother/child`=="m"),]
df.2fields.children <- df.2fields.typings[which(df.2fields.typings$`mother/child`=="c"),]
# Call infer_haplotypes() function using apply() #
df.2fields.haplotypes <- data.frame(matrix(ncol=13, nrow=0))
colnames(df.2fields.haplotypes) <- c("Haplotype", "Sample No.", genes, "Haplotype Length")
df.2fields.haplotypes <- apply(df.2fields.children, 1, infer_haplotypes, df.2fields.mothers, df.2fields.haplotypes)
df.2fields.haplotypes <- do.call(rbind, df.2fields.haplotypes)

# Medium resolution (2 fields plus 'g' group) haplotypes in the whole panel
head(df.4fields.haplotypes)

## Calculate haplotype frequencies and LD metrics in 2 fields + g group inferred haplotypes
BC.2fields.frequencies <- get_frequencies(df.2fields.haplotypes, gene1 = "HLA-B", gene2="HLA-C")
DQB1DRB1.2fields.frequencies <- get_frequencies(df.2fields.haplotypes, gene1 = "HLA-DRB1", gene2="HLA-DQB1")

# BC Frequencies plot
BC.2fields.frequencies <- BC.2fields.frequencies[order(-BC.2fields.frequencies$Haplotype.Frequence),]
BC.2fields.frequencies$rank <- seq_len(nrow(BC.2fields.frequencies))
BC.2fields.frequencies$cumulative <- cumsum(BC.2fields.frequencies$Haplotype.Frequence)

histogram.2fields.BC <- ggplot(BC.2fields.frequencies, aes(x=rank, y=Haplotype.Frequence)) +
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_line(aes(y=cumulative / 10), color="red", size=1.2)+
  scale_y_continuous(name="Haplotype frequency (%)",
                     sec.axis = sec_axis(trans =~.*10, name="Cumulative frequency (%)"))+
  labs(title="HLA-B~C",
       x="Unique haplotype rankings")
histogram.2fields.BC

# DQB1DRB1 Frequences plot
DQB1DRB1.2fields.frequencies <- DQB1DRB1.2fields.frequencies[order(-DQB1DRB1.2fields.frequencies$Haplotype.Frequence),]
DQB1DRB1.2fields.frequencies$rank <- seq_len(nrow(DQB1DRB1.2fields.frequencies))
DQB1DRB1.2fields.frequencies$cumulative <- cumsum(DQB1DRB1.2fields.frequencies$Haplotype.Frequence)

histogram.2fields.DQB1DRB1 <- ggplot(data=DQB1DRB1.2fields.frequencies, aes(x=rank, y=Haplotype.Frequence)) +
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_line(aes(y=cumulative / 10), color="red", size=1.2)+
  scale_y_continuous(name="Haplotype frequency (%)",
                     sec.axis = sec_axis(trans =~.*10, name="Cumulative frequency (%)"))+
  labs(title="HLA-DQB1~DRB1",
       x = "Unique haplotype rankings")
histogram.2fields.DQB1DRB1

# Combine both histograms
grid.arrange(histogram.2fields.BC, histogram.2fields.DQB1DRB1, 
             nrow = 1,  
             top = textGrob("Distributions of two-loci haplotype frequencies in the database", gp=gpar(fontsize=16), vjust = 0.5))

# Plot differences in number of unique haplotypes at 2 fields vs 4 fields resolution
barplot.unique.haplotypes <- data.frame(Haplotype=rep(c("DRB1~DQB1", "B~C"), each=2),
                         Resolution=rep(c("4 fields", "2 fields + 'g' group"),2),
                         Unique_haplotypes=c(nrow(DQB1DRB1.4fields.frequencies),nrow(DQB1DRB1.2fields.frequencies),
                                             nrow(BC.4fields.frequencies),nrow(BC.2fields.frequencies)))
barplot.unique.haplotypes <- ggplot(data=barplot.unique.haplotypes, aes(x=Resolution, y=Unique_haplotypes, fill=Haplotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(
    aes(x = Resolution, y = Unique_haplotypes, label = Unique_haplotypes),
    position = position_dodge(width = 1),
    vjust = -0.4, size = 3)+  
  labs(y="Number of unique haplotypes")+
  theme_bw()
barplot.unique.haplotypes

## Analysis of discordant samples plots
# Numbers were obtained by manual inspection
pie.discordant.genes <- data.frame(
  group = c("DPA1", "G", "DPB1", "DRB1", "DQB1", "DRB345", "A"),
  value = c(8, 1,  10, 2, 3, 3, 1)
)
pie.discordant.genes <- ggplot(pie.discordant.genes, aes(x="", y=value, fill=reorder(group, -value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))+
  labs(fill = "Per HLA gene",
       x = NULL,
       y = NULL,
       title = NULL)+
  theme_minimal()
pie.discordant.genes

pie.discordant.techniques <- data.frame(
  group = c("IMGT version issue", "Allele dropout", "Phasing ambiguity", "Sequence ambiguity"),
  value = c(16, 6, 3, 3)
)
pie.discordant.techniques <- ggplot(pie.discordant.techniques, aes(x="", y=value, fill=reorder(group, -value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))+
  labs(fill = "Per type of error",
       x = NULL,
       y = NULL,
       title = NULL)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
pie.discordant.techniques

# Combine plots showing in which genes and which type of errors led to discordant haplotypes
grid.arrange(pie.discordant.genes, pie.discordant.techniques, 
             nrow = 1,  
             widths = c(0.85,1),
             top = textGrob("Number of typing errors detected by genotype inconsistencies", gp=gpar(fontsize=16), vjust = 3))
