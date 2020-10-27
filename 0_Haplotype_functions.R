####################################################
# Custom functions created for haplotype analysis  #
####################################################

# These libraries are needed for these functions and the rest of the scripts
# Load libraries

library(readxl)
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(heatmaply)
library(ggplot2)
library(gridExtra)
library(grid)
library(lemon)

## Infer haplotypes from related children and mothers HLA genotypes ##
infer_haplotypes <- function(df.children, df.mothers, output){
  
  sample_no <- df.children[["no."]]
  shared_haplotype <- data.frame(haplotype="Shared", sample=sample_no)
  father_haplotype <- data.frame(haplotype="Father", sample=sample_no)
  mother_haplotype <- data.frame(haplotype="Mother", sample=sample_no)
  
  for (gene in genes){
    
    gene2 <- paste(gene, "2") 
    c_alleles <- c(df.children[[gene]], df.children[[gene2]])
    m_alleles <- c(df.mothers[df.mothers$no.==sample_no, gene], df.mothers[df.mothers$no.==sample_no, gene2])
    m_alleles <- unlist(m_alleles, use.names = FALSE)
    shared_alleles <- intersect(c_alleles, m_alleles)
    unshared_father <- setdiff(c_alleles, m_alleles)
    unshared_mother <- setdiff(m_alleles, c_alleles)
    
    if(all(c_alleles=="NA") | all(m_alleles=="NA")){
      if (setequal(c_alleles, m_alleles)){
        father_haplotype[[gene]] <- "NA"
        shared_haplotype[[gene]] <- "NA"
        mother_haplotype[[gene]] <- "NA"
      }
      else{
        father_haplotype[[gene]] <- "Unknown"
        shared_haplotype[[gene]] <- "Unknown"
        mother_haplotype[[gene]] <- "Unknown"
      }    
    }
    else if(setequal(c_alleles, m_alleles)){
      if(all(c(length(unique(c_alleles)), length(unique(m_alleles)))==c(1,1))){
        father_haplotype[[gene]] <- c_alleles[1]
        shared_haplotype[[gene]] <- c_alleles[1]
        mother_haplotype[[gene]] <- c_alleles[1]
      }
      else{
        father_haplotype[[gene]] <- "Ambiguous"
        shared_haplotype[[gene]] <- "Ambiguous"
        mother_haplotype[[gene]] <- "Ambiguous"
      }
    } 
    else{
      if(length(shared_alleles)==0){
        shared_haplotype[[gene]] <- "Discordant"
        father_haplotype[[gene]] <- "Discordant"
        mother_haplotype[[gene]] <- "Discordant"
      }
      else{
        if(length(shared_alleles)==1 & length(unshared_father)==1 & length(unshared_mother)==1){
          shared_haplotype[[gene]] <- shared_alleles
          father_haplotype[[gene]] <- unshared_father
          mother_haplotype[[gene]] <- unshared_mother
        }
        else if(length(shared_alleles)==1 & length(unshared_father)==0 & length(unshared_mother)==1){
          shared_haplotype[[gene]] <- shared_alleles
          father_haplotype[[gene]] <- shared_alleles
          mother_haplotype[[gene]] <- unshared_mother
        }
        else if(length(shared_alleles)==1 & length(unshared_father)==1 & length(unshared_mother)==0){
          shared_haplotype[[gene]] <- shared_alleles
          father_haplotype[[gene]] <- unshared_father
          mother_haplotype[[gene]] <- shared_alleles
        }
        else{
            print("Unknown error")
        }  
      }
    }
  } 

  no_haplotype <- c("NA", "Unknown", "Discordant", "Ambiguous")
  typed_shared <- unlist(shared_haplotype[genes], use.names = FALSE)
  shared_haplotype[["Haplotype Length"]] <- length(typed_shared[! typed_shared %in% no_haplotype])
  typed_mother <- unlist(mother_haplotype[genes], use.names = FALSE)
  mother_haplotype[["Haplotype Length"]] <- length(typed_mother[! typed_mother %in% no_haplotype])
  typed_father <- unlist(father_haplotype[genes], use.names = FALSE)
  father_haplotype[["Haplotype Length"]] <- length(typed_father[! typed_father %in% no_haplotype])
  
  output[nrow(output)+1,] <- shared_haplotype
  output[nrow(output)+1,] <- mother_haplotype
  output[nrow(output)+1,] <- father_haplotype
  return(output)
}

## Calculate 2 loci haplotypic, allelic frequencies ##
## and LD metrics based on a haplotypes dataframe   ##

get_frequencies <- function(df.haplotypes, gene1, gene2){
  
  ##Filter informative haplotypes##
  combined.haplotypes <- (paste0(df.haplotypes[,gene1],",",df.haplotypes[,gene2]))
  combined.haplotypes <-  combined.haplotypes[!grepl("Unknown", combined.haplotypes)]
  combined.haplotypes <-  combined.haplotypes[!grepl("NA", combined.haplotypes)]
  combined.haplotypes <-  combined.haplotypes[!grepl("Discordant", combined.haplotypes)]
  
  ##Get REALLY informative allelic frequencies##
  df.to.alleles <- as.data.frame(combined.haplotypes)
  df.to.alleles <- df.to.alleles %>%
    separate("combined.haplotypes", c(gene1, gene2), ",")
  gene1.allelic.freqs <- df.to.alleles[,gene1]
  gene2.allelic.freqs <- df.to.alleles[,gene2]
  
  ##Calculate allelic frequencies##
  gene1.allelic.freqs <- as.data.frame(sort(table(gene1.allelic.freqs), decreasing = TRUE))
  gene2.allelic.freqs <- as.data.frame(sort(table(gene2.allelic.freqs), decreasing = TRUE))
  
  gene1.allelic.freqs <- gene1.allelic.freqs %>%  
    mutate(Allele1.freq = Freq / sum(Freq))
  gene2.allelic.freqs <- gene2.allelic.freqs %>%  
    mutate(Allele2.freq = Freq / sum(Freq))
  
  colnames(gene1.allelic.freqs)[c(1,2)] <- c(gene1, "Allele1.Count")
  colnames(gene2.allelic.freqs)[c(1,2)] <- c(gene2, "Allele2.Count")
  
  ##Calculate haplotypic frequencies##
  haplotypes.frequencies <- as.data.frame(sort(table(combined.haplotypes), decreasing = TRUE))
  splited.haplotypes <- str_split_fixed(haplotypes.frequencies$combined.haplotypes, ",", 2)
  haplotypes.frequencies <- cbind(splited.haplotypes, haplotypes.frequencies)
  haplotypes.frequencies$combined.haplotypes <- NULL
  colnames(haplotypes.frequencies)[1] <- gene1
  colnames(haplotypes.frequencies)[2] <- gene2
  
  names(haplotypes.frequencies)[names(haplotypes.frequencies) == 'Freq'] <- 'Count'
  haplotypes.frequencies <- haplotypes.frequencies %>%  
    mutate(Haplotype.Frequence = Count / sum(Count))
  
  ##Merge haplotypic and allelic frequencies
  haplotypes.frequencies <- merge(haplotypes.frequencies, gene2.allelic.freqs, by = gene2)
  haplotypes.frequencies <- merge(haplotypes.frequencies, gene1.allelic.freqs, by = gene1)
  
  ##Calculate linkage disequilibrium coefficent (Ddelta)##
  expected <- haplotypes.frequencies[, "Allele1.freq"] * haplotypes.frequencies[, "Allele2.freq"]
  LD <- haplotypes.frequencies$Haplotype.Frequence - expected
  haplotypes.frequencies$D <- LD 
  
  ##Calculate normalized LD, D'##
  normalized.LD <- vector()
  normalized.LD <- apply(haplotypes.frequencies, 1, calculate_normalized_LD, normalized.LD)
  haplotypes.frequencies[["D'"]] <- normalized.LD
  
  ##Calculate allelic proportions in haplotypes##
  haplotypes.frequencies <- haplotypes.frequencies %>%  
    mutate(Allele1.Prop = round(Count / Allele1.Count * 100, digits = 4))
  haplotypes.frequencies <- haplotypes.frequencies %>%  
    mutate(Allele2.Prop = round(Count / Allele2.Count * 100, digits = 4))

  
  ##Normalize and round output frequencies##
  haplotypes.frequencies[, c(4,6,8,9,10)] <- round(haplotypes.frequencies[, c(4,6,8,9,10)] * 100, digits = 2)
  haplotypes.frequencies[, c(11,12)] <- round(haplotypes.frequencies[, c(11,12)], digits = 2 )
  haplotypes.frequencies <- haplotypes.frequencies[c(1,2,3,4,7,8,5,6,9,10,11,12)]
  return(haplotypes.frequencies)
}

## Calculate normalized LD metric (D')
calculate_normalized_LD <- function(haplotypes.frequencies, normalized.LD){
  if (as.numeric(haplotypes.frequencies[["D"]]) >= 0){
    ##Calculate D' for positive values of D##
    dmax <- c((1-as.numeric(haplotypes.frequencies[["Allele2.freq"]]))*as.numeric(haplotypes.frequencies[["Allele1.freq"]]), 
              (1-as.numeric(haplotypes.frequencies[["Allele1.freq"]]))*as.numeric(haplotypes.frequencies[["Allele2.freq"]]))
    dprime <- as.numeric(haplotypes.frequencies[["D"]]) / min(dmax)
  }    
  else{
    ##Calculate D' for negative values of D##
    dmax<- c(-(as.numeric(haplotypes.frequencies[["Allele1.freq"]]) * as.numeric(haplotypes.frequencies[["Allele2.freq"]])),
             -(1-as.numeric(haplotypes.frequencies[["Allele1.freq"]]))*(1-as.numeric(haplotypes.frequencies[["Allele2.freq"]])))
    dprime <- as.numeric(haplotypes.frequencies[["D"]]) / max(dmax)
  }
  normalized.LD <- c(normalized.LD, dprime)
  return(normalized.LD)
}

## Calculate LD metric (D') from special NMDP haplotype frequency table format
calculate_normalized_LD_NMDP <- function(haplotypes.frequencies, normalized.LD, freq1, freq2){
  if (as.numeric(haplotypes.frequencies[["LD"]]) >= 0){
    ##Calculate D' for positive values of D##
    dmax <- c((1-as.numeric(haplotypes.frequencies[[freq2]]))*as.numeric(haplotypes.frequencies[[freq1]]), 
              (1-as.numeric(haplotypes.frequencies[[freq1]]))*as.numeric(haplotypes.frequencies[[freq2]]))
    dprime <- as.numeric(haplotypes.frequencies[["LD"]]) / min(dmax)
    if (is.nan(dprime)){
      dprime <- 0
    }
  }    
  else{
    ##Calculate D' for negative values of D##
    dmax<- c(-(as.numeric(haplotypes.frequencies[[freq1]]) * as.numeric(haplotypes.frequencies[[freq2]])),
             -(1-as.numeric(haplotypes.frequencies[[freq1]]))*(1-as.numeric(haplotypes.frequencies[[freq2]])))
    dprime <- as.numeric(haplotypes.frequencies[["LD"]]) / max(dmax)
  }
  normalized.LD <- c(normalized.LD, dprime)
  return(normalized.LD)
}

## Translate two fields typings to two fields plus 'g' group nomenclature
translate_g_group <- function(df.typings, translation){
  df.typings <- paste0(df.typings, ",")
  translation$Alleles <- paste0(translation$Alleles, ",")
  if (TRUE %in% (grepl(df.typings, translation$Alleles, fixed=TRUE))){
    df.typings <- translation[grepl(df.typings, translation$Alleles, fixed=TRUE),]$Gallele
  }
  else if(grepl("NA", df.typings, fixed=TRUE)){
    df.typings <- "NA"
  }
  else{
    df.typings <- str_replace(df.typings, ",", "")
  }
  return(df.typings)
}

## Reduce typing resolution from four fields to two fields
reduce_typing_resolution <- function(df.haplotypes, resolution){
  if (resolution==1){
    apply(df.haplotypes, 2, function(df.haplotypes){
      new.df.haplotypes <- str_split(df.haplotypes, ":")
      lapply(new.df.haplotypes, function(new.df.haplotypes){
        if (length(new.df.haplotypes) > resolution){
          new.df.haplotypes <- new.df.haplotypes[c(1)]
        }
        else{
          new.df.haplotypes <- new.df.haplotypes
        }
      })
    })
  }
  else {
    apply(df.haplotypes, 2, function(df.haplotypes){
      new.df.haplotypes <- str_split(df.haplotypes, ":")
      lapply(new.df.haplotypes, function(new.df.haplotypes){
        if (length(new.df.haplotypes) >= resolution){
          #new.df.haplotypes <- new.df.haplotypes[c(1,2)]
          new.df.haplotypes <- paste0(new.df.haplotypes[1], ":", new.df.haplotypes[2])
        }
        else{
          new.df.haplotypes <- new.df.haplotypes
        }
      })
    })
  }
}

## Estimate combinations of two loci haplotypes in a single sample typing, using
## NMDP frequency tables. Calculate estimation metrics
estimate_haplotypes <- function(sample.typing, gene1, gene2, freq.table, output){
  gene1.alleles <- c(sample.typing[gene1], sample.typing[paste0(gene1, " 2")])
  gene2.alleles <- c(sample.typing[gene2], sample.typing[paste0(gene2, " 2")])
  all.alleles <- c(gene1.alleles, gene2.alleles)
  if (length(unique(all.alleles))== 2){
    possible.haplotypes <- rbind(crossing(gene1.alleles, gene2.alleles), crossing(gene1.alleles, gene2.alleles))
    possible.haplotypes <- possible.haplotypes %>% unite(combined, 1, 2, remove = FALSE)
  }
  else{
  possible.haplotypes <- (crossing(gene1.alleles, gene2.alleles))
  possible.haplotypes <- possible.haplotypes %>% unite(combined, 1, 2, remove = FALSE)
  }
  possible.haplotypes <- merge(possible.haplotypes, freq.table, by="combined", all.x=TRUE)
  comb.list <- combn(possible.haplotypes$combined, 2)
  comb.list <- as.data.frame(t(comb.list))
  comb.list$alleles <- paste(comb.list[,1], comb.list[,2], sep="_")
  comb.list$alleles <- str_split(comb.list$alleles, "_")
  
  comb.list <- comb.list %>% 
    mutate(real = lapply(comb.list$alleles, setequal, all.alleles))
  comb.list <- comb.list[which(comb.list$real==TRUE),]
  comb.list <- merge(comb.list, possible.haplotypes, by.x = "V1", by.y = "combined")
  comb.list <- merge(comb.list, possible.haplotypes, by.x = "V2", by.y= "combined")
  comb.list <- comb.list %>%
    mutate(metric = CAU_freq.x * CAU_freq.y)
  comb.list$normalized <- (comb.list$metric / max(comb.list$metric))
  comb.list$sample <- sample.typing["no."]
  comb.list$`mother/child` <- sample.typing["mother/child"]
  comb.list$unique.alleles <- length(unique(all.alleles))
  if (length(unique(all.alleles))== 2){
    comb.list <- comb.list[1,]
  }
  output <- rbind(output, comb.list)

  return(output)
  }

## Allele exchanging method for homozygous genotypes
allele_exchanging_homozygous <- function(likely_haplotype, freq.table, output){
  likely_haplotype <- as.list(likely_haplotype)
  likely_haplotype$CAU_rank.x <- as.numeric(likely_haplotype$CAU_rank.x)
  likely_haplotype$CAU_rank.y <- as.numeric(likely_haplotype$CAU_rank.y)
  
  if(likely_haplotype$gene1.alleles.x == likely_haplotype$gene1.alleles.y){
    allele2 <- freq.table[grepl(likely_haplotype$gene2.alleles.x, freq.table$combined, fixed=TRUE),]
    allele2$prev.haplotype <- likely_haplotype$haplotype.x
    allele2$prev.allele <- likely_haplotype$gene1.alleles.x
    allele2$prev.CAU_freq <- likely_haplotype$CAU_freq.x
    allele2$prev.CAU_rank <- likely_haplotype$CAU_rank.x
    allele2$prev.LD <- likely_haplotype$LD.x
    allele2$`prev.D'` <- likely_haplotype$`D'.x`
    allele2$prev.allele1prop <- likely_haplotype$allele1prop.x
    allele2$prev.allele2prop <- likely_haplotype$allele2prop.x
    allele2$sample <- likely_haplotype$sample
    allele2$`mother/child` <- likely_haplotype$`mother/child`
    allele2$unique.alleles <- likely_haplotype$unique.alleles
    allele2 <- allele2[which(allele2$CAU_rank < likely_haplotype$CAU_rank.x),]
    
    top3 <- tail(sort(allele2$allele1prop), 3)
    allele2 <- allele2[which(allele2$allele1prop %in% top3),]
    if(nrow(allele2)>0){
    allele2$propsum <- sum(allele2$allele1prop)
    }
    
    allele4 <- freq.table[grepl(likely_haplotype$gene2.alleles.y, freq.table$combined, fixed = TRUE),]
    allele4$prev.allele <- likely_haplotype$gene1.alleles.y
    allele4$prev.haplotype <- likely_haplotype$haplotype.y
    allele4$prev.CAU_freq <- likely_haplotype$CAU_freq.y
    allele4$prev.CAU_rank <- likely_haplotype$CAU_rank.y
    allele4$prev.LD <- likely_haplotype$LD.y
    allele4$`prev.D'` <- likely_haplotype$`D'.y`
    allele4$prev.allele1prop <- likely_haplotype$allele1prop.y
    allele4$prev.allele2prop <- likely_haplotype$allele2prop.y
    allele4$sample <- likely_haplotype$sample
    allele4$`mother/child` <- likely_haplotype$`mother/child`
    allele4$unique.alleles <- likely_haplotype$unique.alleles
    allele4 <- allele4[which(allele4$CAU_rank < likely_haplotype$CAU_rank.y),]
    
    top3 <- tail(sort(allele4$allele1prop),3)
    allele4 <- allele4[which(allele4$allele1prop %in%top3),]
    if(nrow(allele4)>0){
    allele4$propsum <- sum(allele4$allele1prop)
    }
    
    sample.output  <- rbind(allele2, allele4)
    output <- rbind(output, sample.output)
  }
  
  if(likely_haplotype$gene2.alleles.x == likely_haplotype$gene2.alleles.y){
    allele1 <- freq.table[grepl(likely_haplotype$gene1.alleles.x, freq.table$combined, fixed = TRUE),]
    allele1$prev.haplotype <- likely_haplotype$haplotype.x
    allele1$prev.allele <- likely_haplotype$gene2.alleles.x
    allele1$prev.CAU_freq <- likely_haplotype$CAU_freq.x
    allele1$prev.CAU_rank <- likely_haplotype$CAU_rank.x
    allele1$prev.LD <- likely_haplotype$LD.x
    allele1$`prev.D'` <- likely_haplotype$`D'.x`
    allele1$prev.allele1prop <- likely_haplotype$allele1prop.x
    allele1$prev.allele2prop <- likely_haplotype$allele2prop.x
    allele1$sample <- likely_haplotype$sample
    allele1$`mother/child` <- likely_haplotype$`mother/child`
    allele1$unique.alleles <- likely_haplotype$unique.alleles
    allele1 <- allele1[which(allele1$CAU_rank < likely_haplotype$CAU_rank.x),]
    
    top3 <- tail(sort(allele1$allele2prop),3)
    allele1 <- allele1[which(allele1$allele2prop %in%top3),]
    if(nrow(allele1)>0){
    allele1$propsum <- sum(allele1$allele2prop)
    }
    
    allele3 <- freq.table[grepl(likely_haplotype$gene1.alleles.y, freq.table$combined, fixed=TRUE),]
    allele3$prev.allele <- likely_haplotype$gene2.alleles.y
    allele3$prev.haplotype <- likely_haplotype$haplotype.y
    allele3$prev.CAU_freq <- likely_haplotype$CAU_freq.y
    allele3$prev.CAU_rank <- likely_haplotype$CAU_rank.y
    allele3$prev.LD <- likely_haplotype$LD.y
    allele3$`prev.D'` <- likely_haplotype$`D'.y`
    allele3$prev.allele1prop <- likely_haplotype$allele1prop.y
    allele3$prev.allele2prop <- likely_haplotype$allele2prop.y
    allele3$sample <- likely_haplotype$sample
    allele3$`mother/child` <- likely_haplotype$`mother/child`
    allele3$unique.alleles <- likely_haplotype$unique.alleles
    allele3 <- allele3[which(allele3$CAU_rank < likely_haplotype$CAU_rank.y),]

    top3 <- tail(sort(allele3$allele2prop), 3)
    allele3 <- allele3[which(allele3$allele2prop %in% top3),]
    if(nrow(allele3)>0){
    allele3$propsum <- sum(allele3$allele2prop)
    }
    
    sample.output  <- rbind(allele1, allele3)
    output <- rbind(output, sample.output)
  }
  return(output)
}

##Allele exchanging method for complete alleles
allele_exchanging_haplotypes <- function(likely_haplotype, freq.table, output) {
  likely_haplotype <- as.list(likely_haplotype)
  likely_haplotype$CAU_rank.x <- as.numeric(likely_haplotype$CAU_rank.x)
  likely_haplotype$CAU_rank.y <- as.numeric(likely_haplotype$CAU_rank.y)
  
  allele1 <- freq.table[grepl(likely_haplotype$gene1.alleles.x, freq.table$combined, fixed = TRUE),]
  allele1$prev.haplotype <- likely_haplotype$haplotype.x
  allele1$prev.allele <- likely_haplotype$gene2.alleles.x
  allele1$prev.CAU_freq <- likely_haplotype$CAU_freq.x
  allele1$prev.CAU_rank <- likely_haplotype$CAU_rank.x
  allele1$prev.LD <- likely_haplotype$LD.x
  allele1$`prev.D'` <- likely_haplotype$`D'.x`
  allele1$prev.allele1prop <- likely_haplotype$allele1prop.x
  allele1$prev.allele2prop <- likely_haplotype$allele2prop.x
  allele1$sample <- likely_haplotype$sample
  allele1$`mother/child` <- likely_haplotype$`mother/child`
  allele1$unique.alleles <- likely_haplotype$unique.alleles
  allele1 <- allele1[which(allele1$CAU_rank < likely_haplotype$CAU_rank.x),]

  allele2 <- freq.table[grepl(likely_haplotype$gene2.alleles.x, freq.table$combined, fixed=TRUE),]
  allele2$prev.haplotype <- likely_haplotype$haplotype.x
  allele2$prev.allele <- likely_haplotype$gene1.alleles.x
  allele2$prev.CAU_freq <- likely_haplotype$CAU_freq.x
  allele2$prev.CAU_rank <- likely_haplotype$CAU_rank.x
  allele2$prev.LD <- likely_haplotype$LD.x
  allele2$`prev.D'` <- likely_haplotype$`D'.x`
  allele2$prev.allele1prop <- likely_haplotype$allele1prop.x
  allele2$prev.allele2prop <- likely_haplotype$allele2prop.x
  allele2$sample <- likely_haplotype$sample
  allele2$`mother/child` <- likely_haplotype$`mother/child`
  allele2$unique.alleles <- likely_haplotype$unique.alleles
  allele2 <- allele2[which(allele2$CAU_rank < likely_haplotype$CAU_rank.x),]

  allele3 <- freq.table[grepl(likely_haplotype$gene1.alleles.y, freq.table$combined, fixed=TRUE),]
  allele3$prev.allele <- likely_haplotype$gene2.alleles.y
  allele3$prev.haplotype <- likely_haplotype$haplotype.y
  allele3$prev.CAU_freq <- likely_haplotype$CAU_freq.y
  allele3$prev.CAU_rank <- likely_haplotype$CAU_rank.y
  allele3$prev.LD <- likely_haplotype$LD.y
  allele3$`prev.D'` <- likely_haplotype$`D'.y`
  allele3$prev.allele1prop <- likely_haplotype$allele1prop.y
  allele3$prev.allele2prop <- likely_haplotype$allele2prop.y
  allele3$sample <- likely_haplotype$sample
  allele3$`mother/child` <- likely_haplotype$`mother/child`
  allele3$unique.alleles <- likely_haplotype$unique.alleles
  allele3 <- allele3[which(allele3$CAU_rank < likely_haplotype$CAU_rank.y),]

  allele4 <- freq.table[grepl(likely_haplotype$gene2.alleles.y, freq.table$combined, fixed = TRUE),]
  allele4$prev.allele <- likely_haplotype$gene1.alleles.y
  allele4$prev.haplotype <- likely_haplotype$haplotype.y
  allele4$prev.CAU_freq <- likely_haplotype$CAU_freq.y
  allele4$prev.CAU_rank <- likely_haplotype$CAU_rank.y
  allele4$prev.LD <- likely_haplotype$LD.y
  allele4$`prev.D'` <- likely_haplotype$`D'.y`
  allele4$prev.allele1prop <- likely_haplotype$allele1prop.y
  allele4$prev.allele2prop <- likely_haplotype$allele2prop.y
  allele4$sample <- likely_haplotype$sample
  allele4$`mother/child` <- likely_haplotype$`mother/child`
  allele4$unique.alleles <- likely_haplotype$unique.alleles
  allele4 <- allele4[which(allele4$CAU_rank < likely_haplotype$CAU_rank.y),]

  output  <- rbind(allele1, allele2, allele3, allele4)
  
  return(output)
}

##Allele exchanging method for allele subtypes
allele_exchanging_subtypes <- function(likely_haplotype, freq.table, output) {
  likely_haplotype <- as.list(likely_haplotype)
  likely_haplotype$CAU_rank.x <- as.numeric(likely_haplotype$CAU_rank.x)
  likely_haplotype$CAU_rank.y <- as.numeric(likely_haplotype$CAU_rank.y)

  subtype1 <- paste(likely_haplotype$gene1.alleles.x, 
                    unlist(str_split(likely_haplotype$gene2.alleles.x, ":"))[1], sep="_")
  subtype2 <- paste(likely_haplotype$gene2.alleles.x, 
                    unlist(str_split(likely_haplotype$gene1.alleles.x, ":"))[1], sep="_")
  
  subtype3 <- paste(likely_haplotype$gene1.alleles.y, 
                    unlist(str_split(likely_haplotype$gene2.alleles.y, ":"))[1], sep="_")
  subtype4 <- paste(likely_haplotype$gene2.alleles.y, 
                    unlist(str_split(likely_haplotype$gene1.alleles.y, ":"))[1], sep="_")

  allele1 <- freq.table[grepl(subtype1, freq.table$combined, fixed=TRUE),]
  allele1$prev.haplotype <- likely_haplotype$haplotype.x
  allele1$prev.allele <- likely_haplotype$gene2.alleles.x
  allele1$prev.CAU_freq <- likely_haplotype$CAU_freq.x
  allele1$prev.CAU_rank <- likely_haplotype$CAU_rank.x
  allele1$prev.LD <- likely_haplotype$LD.x
  allele1$`prev.D'` <- likely_haplotype$`D'.x`
  allele1$prev.allele1prop <- likely_haplotype$allele1prop.x
  allele1$prev.allele2prop <- likely_haplotype$allele2prop.x
  allele1$sample <- likely_haplotype$sample
  allele1$`mother/child` <- likely_haplotype$`mother/child`
  allele1$unique.alleles <- likely_haplotype$unique.alleles
  allele1 <- allele1[which(allele1$CAU_rank < likely_haplotype$CAU_rank.x),]

  allele2 <- freq.table[grepl(subtype2, freq.table$combined_inverse, fixed=TRUE),]
  allele2$prev.haplotype <- likely_haplotype$haplotype.x
  allele2$prev.allele <- likely_haplotype$gene1.alleles.x
  allele2$prev.CAU_freq <- likely_haplotype$CAU_freq.x
  allele2$prev.CAU_rank <- likely_haplotype$CAU_rank.x
  allele2$prev.LD <- likely_haplotype$LD.x
  allele2$`prev.D'` <- likely_haplotype$`D'.x`
  allele2$prev.allele1prop <- likely_haplotype$allele1prop.x
  allele2$prev.allele2prop <- likely_haplotype$allele2prop.x
  allele2$sample <- likely_haplotype$sample
  allele2$`mother/child` <- likely_haplotype$`mother/child`
  allele2$unique.alleles <- likely_haplotype$unique.alleles
  allele2 <- allele2[which(allele2$CAU_rank < likely_haplotype$CAU_rank.x),]

  allele3 <- freq.table[grepl(subtype3, freq.table$combined, fixed=TRUE),]
  allele3$prev.allele <- likely_haplotype$gene2.alleles.y
  allele3$prev.haplotype <- likely_haplotype$haplotype.y
  allele3$prev.CAU_freq <- likely_haplotype$CAU_freq.y
  allele3$prev.CAU_rank <- likely_haplotype$CAU_rank.y
  allele3$prev.LD <- likely_haplotype$LD.y
  allele3$`prev.D'` <- likely_haplotype$`D'.y`
  allele3$prev.allele1prop <- likely_haplotype$allele1prop.y
  allele3$prev.allele2prop <- likely_haplotype$allele2prop.y
  allele3$sample <- likely_haplotype$sample
  allele3$`mother/child` <- likely_haplotype$`mother/child`
  allele3$unique.alleles <- likely_haplotype$unique.alleles
  allele3 <- allele3[which(allele3$CAU_rank < likely_haplotype$CAU_rank.y),]

  allele4 <- freq.table[grepl(subtype4, freq.table$combined_inverse, fixed = TRUE),]
  allele4$prev.allele <- likely_haplotype$gene1.alleles.y
  allele4$prev.haplotype <- likely_haplotype$haplotype.y
  allele4$prev.CAU_freq <- likely_haplotype$CAU_freq.y
  allele4$prev.CAU_rank <- likely_haplotype$CAU_rank.y
  allele4$prev.LD <- likely_haplotype$LD.y
  allele4$`prev.D'` <- likely_haplotype$`D'.y`
  allele4$prev.allele1prop <- likely_haplotype$allele1prop.y
  allele4$prev.allele2prop <- likely_haplotype$allele2prop.y
  allele4$sample <- likely_haplotype$sample
  allele4$`mother/child` <- likely_haplotype$`mother/child`
  allele4$unique.alleles <- likely_haplotype$unique.alleles
  allele4 <- allele4[which(allele4$CAU_rank < likely_haplotype$CAU_rank.y),]

  output <- rbind(allele1, allele2, allele3, allele4)
  
  return(output)
}

## Validation of estimated haplotypes, by comparison with inferred haplotypes
validate_haplotypes <- function(observed.haplotypes, estimated.haplotypes, output.val){
  output.val <- data.frame(matrix(ncol=4, nrow=0))
  output.val <- apply(estimated.haplotypes, 1, validate, observed.haplotypes, output.val)
  output.val <- do.call(rbind, output.val)
  
  return(output.val)
}
# Subordinated function inside validate_haplotypes()
validate <- function(estimated.haplotypes, observed.haplotypes, output.val){
  sample <- estimated.haplotypes["sample"]
  observed <- observed.haplotypes[which(observed.haplotypes$sample==estimated.haplotypes["sample"]), "combined"]
  validation.x <- list()
  validation.y <- list()
  validation.x["sample"] <- sample
  validation.y["sample"] <- sample 
  validation.x["haplotype"] <- estimated.haplotypes["haplotype.x"]
  validation.y["haplotype"] <- estimated.haplotypes["haplotype.y"]
  validation.x["observed"] <- paste0(observed[1], "+", observed[2])
  validation.y["observed"] <- paste0(observed[1], "+", observed[2])
  if(grepl("Ambiguous", observed) | grepl("Unknown", observed) | grepl("Discordant", observed)){
    validation.x["validated"] <- "Not validated"
    validation.y["validated"] <- "Not validated"
  }
  else{
    validation.x["validated"] <- estimated.haplotypes["haplotype.x"] %in% observed
    validation.y["validated"] <- estimated.haplotypes["haplotype.y"] %in% observed 
  }
  output.val <- rbind(output.val, validation.x, validation.y)
  
  return(output.val)
}
