##################################################################
## Correct typing results that were re analyzed after discovery ##
## of typing error upon inconsistencies in the haplotypes       ##
##################################################################

## If these typing results are corrected, different number of inferred haplotypes,
# estimation results, allele exchanging method, will be obtained than those reported in the
# research project

df.sample.overview[which(df.sample.overview$`no.`=="26" &
                           df.sample.overview$`mother/child`=="c"), "HLA-DPB1 2"] <- "124:01:01:01"
df.sample.overview[which(df.sample.overview$`no.`=="26" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPB1"] <- "4:01:01:03"

df.sample.overview[which(df.sample.overview$`no.`=="31" &
                           df.sample.overview$`mother/child`=="c"), "HLA-DPA1 2"] <- "01:03:01:01"

df.sample.overview[which(df.sample.overview$`no.`=="44" &
                           df.sample.overview$`mother/child`=="c"), "HLA-DPA1"] <- "01:03:01:05"

df.sample.overview[which(df.sample.overview$`no.`=="48" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DQB1"] <- "03:01:01:02"

df.sample.overview[which(df.sample.overview$no.=="61"), c("HLA-G", "HLA-G 2")] <- "NA"

df.sample.overview[which(df.sample.overview$`no.`=="66" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DRB345"] <- "4*01:03:01:01"

df.sample.overview[which(df.sample.overview$no.=="105"), c("HLA-DQB1", "HLA-DQB1 2")] <- "NA"

df.sample.overview[which(df.sample.overview$`no.`=="109"), "HLA-DPB1"] <- "1082:01"

df.sample.overview[which(df.sample.overview$`no.`=="130" &
                           df.sample.overview$`mother/child`=="c"),"HLA-DQB1"] <- "03:01:01:02"

df.sample.overview[which(df.sample.overview$`no.`=="134" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPA1 2"] <- "01:03:01:05"

df.sample.overview[which(df.sample.overview$`no.`=="180" &
                           df.sample.overview$`mother/child`=="c"), "HLA-DRB1 2"] <- "14:54:01:01"

df.sample.overview[which(df.sample.overview$`no.`=="200" &
                           df.sample.overview$`mother/child`=="m"), "HLA-A 2"] <- "68:02:01:01"

df.sample.overview[which(df.sample.overview$`no.`=="244" &
                           df.sample.overview$`mother/child`=="c"),"HLA-DPB1"] <- "09:01:01"
df.sample.overview[which(df.sample.overview$`no.`=="244" &
                           df.sample.overview$`mother/child`=="c"),"HLA-DPB1 2"] <- "02:01:02:01"

df.sample.overview[which(df.sample.overview$no.=="292"), c("HLA-DPB1", "HLA-DPB1 2")] <- "NA"

df.sample.overview[which(df.sample.overview$`no.`=="301" &
                           df.sample.overview$`mother/child`=="c"), "HLA-DPB1"] <- "03:01:01:04"
df.sample.overview[which(df.sample.overview$`no.`=="301" &
                           df.sample.overview$`mother/child`=="c"),"HLA-DPB1 2"] <- "13:01:01:02"

df.sample.overview[which(df.sample.overview$`no.`=="305" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPA1"] <- "01:03:01:05"
df.sample.overview[which(df.sample.overview$`no.`=="319" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPA1"] <- "01:03:01:05"
df.sample.overview[which(df.sample.overview$`no.`=="327" &
                           df.sample.overview$`mother/child`=="c"), "HLA-DPA1 2"] <- "01:03:01:05"
df.sample.overview[which(df.sample.overview$`no.`=="346" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPA1"] <- "01:03:01:05"

df.sample.overview[which(df.sample.overview$`no.`=="373" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPB1"] <- "04:01:01:01"
df.sample.overview[which(df.sample.overview$`no.`=="373" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPB1 2"] <- "04:01:01:02" # might be new allele

df.sample.overview[which(df.sample.overview$no.=="339"), "HLA-DPB1"] <- "02:01:02:01"
df.sample.overview[which(df.sample.overview$no.=="339"), "HLA-DPB1 2"] <- "849:01"

df.sample.overview[which(df.sample.overview$`no.`=="435" &
                           df.sample.overview$`mother/child`=="m"), "HLA-DPA1"] <- "01:03:01:02"
