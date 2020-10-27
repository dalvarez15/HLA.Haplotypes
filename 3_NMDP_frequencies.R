##############################################
## HLA-DRB1~DQB1 NMDP haplotype frequencies  #
## and LD metrics                            #
##############################################

# Import HLA-DRB1~DQB1 NMDP haplotype frequencies tables 
nmdp.DQB1DRB1 <- read_excel("Input/DRB1_DQB1.xlsx")
nmdp.DQB1DRB1[nmdp.DQB1DRB1=="NA"] <- NA
nmdp.DQB1DRB1[, seq(4,50, by=2)] <- sapply(nmdp.DQB1DRB1[,seq(4,50, by=2)], as.numeric)
nmdp.DQB1DRB1 <- nmdp.DQB1DRB1 %>%
  select(DRB1, DQB1, CAU_freq, CAU_rank) #select CAU population
nmdp.DQB1DRB1$combined <- paste0(nmdp.DQB1DRB1$DRB1, "_", nmdp.DQB1DRB1$DQB1)
nmdp.DQB1DRB1 <- nmdp.DQB1DRB1 %>%
  mutate(combined_inverse = paste0(nmdp.DQB1DRB1$DQB1, "_", nmdp.DQB1DRB1$DRB1))

# Import HLA-DRB1 and HLA-DQB1 NMDP allele frequencies tables
nmdp.DRB1 <- read_excel("Input/DRB1.xlsx")
nmdp.DRB1[nmdp.DRB1=="NA"] <- NA
nmdp.DRB1[, seq(2,50, by=2)] <- sapply(nmdp.DRB1[,seq(2,50, by=2)], as.numeric)
nmdp.DRB1 <- nmdp.DRB1 %>%
  select(DRB1, CAU_freq, CAU_rank) #select broad race groups

nmdp.DQB1 <- read_excel("Input/DQB1.xlsx")
nmdp.DQB1[nmdp.DQB1=="NA"] <- NA
nmdp.DQB1[, seq(2,50, by=2)] <- sapply(nmdp.DQB1[,seq(2,50, by=2)], as.numeric)
nmdp.DQB1 <- nmdp.DQB1 %>%
  select(DQB1, CAU_freq, CAU_rank) #select broad race groups

# Combine allele and haplotype frequencies
nmdp.DQB1DRB1 <- merge(nmdp.DQB1DRB1, nmdp.DRB1, by = "DRB1")
nmdp.DQB1DRB1 <- merge(nmdp.DQB1DRB1, nmdp.DQB1, by = "DQB1")
colnames(nmdp.DQB1DRB1)[c(3,4,7,8,9,10)] <- c("CAU_freq", "CAU_rank", "DRB1.CAU_freq",
                                              "DRB1.CAU_rank", "DQB1.CAU_freq", "DQB1.CAU_rank")
## Calculate LD metrics 
# Calculate LD coefficients (D and D')
nmdp.DQB1DRB1 <- nmdp.DQB1DRB1 %>%
  mutate(LD = CAU_freq - (DRB1.CAU_freq*DQB1.CAU_freq))
nmdp.DQB1DRB1[["D'"]] <- apply(nmdp.DQB1DRB1, 1, calculate_normalized_LD_NMDP, 
  nmdp.DQB1DRB1[["D'"]] , "DQB1.CAU_freq", "DRB1.CAU_freq")

# Calculate Allele1prop and Allele2prop
nmdp.DQB1DRB1 <- nmdp.DQB1DRB1 %>%
  mutate(allele1prop = CAU_freq / DQB1.CAU_freq *100)
nmdp.DQB1DRB1 <- nmdp.DQB1DRB1 %>%
  mutate(allele2prop = CAU_freq /DRB1.CAU_freq *100)

# HLA-DQB1~DRB1 haplotype frequencies estimated by the NMDP in Caucasian population
# With LD metrics
head(nmdp.DQB1DRB1[order(nmdp.DQB1DRB1$CAU_rank),])

########################################
## HLA-B~C NMDP haplotype frequencies ##
## and LD metrics                     ##
########################################

# Import NMDP HLA-B~C haplotype frequencies tables 
nmdp.BC <- read_excel("Input/C_B.xlsx")
nmdp.BC[nmdp.BC=="NA"] <- NA
nmdp.BC[, seq(4,54, by=2)] <- sapply(nmdp.BC[,seq(4,54, by=2)], as.numeric)
nmdp.BC <- nmdp.BC %>%
  select(C, B, CAU_freq, CAU_rank)
nmdp.BC <- nmdp.BC %>%
  select(C, B, CAU_freq, CAU_rank)  
nmdp.BC$combined <- paste0(nmdp.BC$C, "_", nmdp.BC$B)
nmdp.BC <- nmdp.BC %>%
  mutate(combined_inverse = paste0(nmdp.BC$B, "_", nmdp.BC$C))

# Import HLA-B and HLA-C NMDP allele frequencies
nmdp.B <- read_excel("Input/B.xlsx")
nmdp.B[nmdp.B=="NA"] <- NA
nmdp.B[, seq(2,50, by=2)] <- sapply(nmdp.B[,seq(2,50, by=2)], as.numeric)
nmdp.B<- nmdp.B %>%
  select(B, CAU_freq, CAU_rank) 
nmdp.C <- read_excel("Input/C.xlsx")
nmdp.C[nmdp.C=="NA"] <- NA
nmdp.C[, seq(2,50, by=2)] <- sapply(nmdp.C[,seq(2,50, by=2)], as.numeric)
nmdp.C <- nmdp.C %>%
  select(C, CAU_freq, CAU_rank) 

# Combine allele and haplotype frequencies
nmdp.BC <- merge(nmdp.BC, nmdp.B, by = "B")
nmdp.BC <- merge(nmdp.BC, nmdp.C, by = "C")
colnames(nmdp.BC)[c(3,4,7,8,9,10)] <- c("CAU_freq", "CAU_rank", "C.CAU_freq",
                                        "C.CAU_rank", "B.CAU_freq", "B.CAU_rank")
## Calculate LD metrics 
# Calculate LD coefficients (D and D')
nmdp.BC <- nmdp.BC %>%
  mutate(LD = CAU_freq - (B.CAU_freq*C.CAU_freq))
nmdp.BC[["D'"]] <- apply(nmdp.BC, 1, calculate_normalized_LD_NMDP, nmdp.BC[["D'"]], 
                         "B.CAU_freq", "C.CAU_freq")
# Calculate Allele1prop and Allele2prop
nmdp.BC <- nmdp.BC %>%
  mutate(allele1prop = CAU_freq / C.CAU_freq *100)
nmdp.BC <- nmdp.BC %>%
  mutate(allele2prop = CAU_freq / B.CAU_freq *100)

# HLA-B~C haplotype frequencies estimated by the NMDP in Caucasian population
# With LD metrics
head(nmdp.BC[order(nmdp.BC$CAU_rank),])
