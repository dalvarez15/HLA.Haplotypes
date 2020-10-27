######################################################
# Read excel files containing HLA typing results for #
# Basel's panel with mother-child couples            #
# Curate database into data.frame format             #
######################################################

## Read excel file using library(readxl)
df.sample.overview <- read_excel("Input/NGSSampleOverview_Basel.xlsx")

## Fix input: variables/genes names, rejected and NA values
names(df.sample.overview)[seq(6, 28, by=2)] <- paste(names(df.sample.overview)[seq(6, 28, by=2)-1], "2")
df.sample.overview[df.sample.overview=="N.a."] <- NA
df.sample.overview <- df.sample.overview[-c(15, 16, 17),]
df.sample.overview[df.sample.overview=="Rejected"] <- NA

# Combine HLA-DRB3,4 and 5 typings into a single gene locus
df.sample.overview <- df.sample.overview %>%
  unite("HLA-DRB345", "HLA-DRB3", "HLA-DRB3 2", 
      "HLA-DRB4", "HLA-DRB4 2","HLA-DRB5", "HLA-DRB5 2", na.rm = TRUE)
df.sample.overview <- df.sample.overview %>%
  separate("HLA-DRB345", c("HLA-DRB3452", "HLA-DRB345"), "_")
colnames(df.sample.overview)[c(13,14)] <- c("HLA-DRB345", "HLA-DRB345 2")
df.sample.overview[df.sample.overview==""] <- "NA"

# Correct NA values and typos
df.sample.overview[is.na(df.sample.overview)] <- "NA"
df.sample.overview[which(df.sample.overview$no. == "394->395"),"no."] <- "395"
names(df.sample.overview)[names(df.sample.overview) == "mother/cild"] <- "mother/child"

# Full resolution typings of all samples in df.sample.overwiew
head(df.sample.overview)
