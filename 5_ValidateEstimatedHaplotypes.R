#############################################################
## Validation of estimated 2-loci haplotypes by comparison ##
## with inferred haplotypes                                ##
#############################################################

## Get observed haplotypes from children and mothers (DRB1, DQB1)
df.2fields.haplotypes.children.DQDR <- subset(df.2fields.haplotypes, `Haplotype`=="Shared" | `Haplotype`=="Father") %>%
  select("Sample No.", "HLA-DQB1", "HLA-DRB1")
df.2fields.haplotypes.children.DQDR$combined <- paste0(df.2fields.haplotypes.children.DQDR[,3], "_", df.2fields.haplotypes.children.DQDR[,2])
colnames(df.2fields.haplotypes.children.DQDR)[1] <- "sample"

df.2fields.haplotypes.mothers.DQDR <- subset(df.2fields.haplotypes, `Haplotype`=="Shared" | `Haplotype`=="Mother") %>%
  select("Sample No.", "HLA-DQB1", "HLA-DRB1")
df.2fields.haplotypes.mothers.DQDR$combined <- paste0(df.2fields.haplotypes.mothers.DQDR[,3], "_", df.2fields.haplotypes.mothers.DQDR[,2])
colnames(df.2fields.haplotypes.mothers.DQDR)[1] <- "sample"

## Get observed haplotypes from children and mothers (B, C)
df.2fields.haplotypes.children.BC <- subset(df.2fields.haplotypes, `Haplotype`=="Shared" | `Haplotype`=="Father") %>%
  select("Sample No.", "HLA-C", "HLA-B")
df.2fields.haplotypes.children.BC$combined <- paste0(df.2fields.haplotypes.children.BC[,2], "_", df.2fields.haplotypes.children.BC[,3])
colnames(df.2fields.haplotypes.children.BC)[1] <- "sample"

df.2fields.haplotypes.mothers.BC <- subset(df.2fields.haplotypes, `Haplotype`=="Shared" | `Haplotype`=="Mother") %>%
  select("Sample No.", "HLA-C", "HLA-B")
df.2fields.haplotypes.mothers.BC$combined <- paste0(df.2fields.haplotypes.mothers.BC[,2], "_", df.2fields.haplotypes.mothers.BC[,3])
colnames(df.2fields.haplotypes.mothers.BC)[1] <- "sample"

## Validate  using validate_haplotypes() function 
validation.BC.c <- validate_haplotypes(df.2fields.haplotypes.children.BC, subset(BC.estimated.haplotypes, `mother/child`=="c"), validation.BC.c)
validation.BC.m <- validate_haplotypes(df.2fields.haplotypes.mothers.BC, subset(BC.estimated.haplotypes, `mother/child`=="m"), validation.BC.m)
validation.DQB1DRB1.m <- validate_haplotypes(df.2fields.haplotypes.mothers.DQDR, subset(DQDR.estimated.haplotypes, `mother/child`=="m"), validation.DQB1DRB1.m)
validation.DQB1DRB1.c <- validate_haplotypes(df.2fields.haplotypes.children.DQDR, subset(DQDR.estimated.haplotypes, `mother/child`=="c"), validation.DQB1DRB1.c)

# Convert validation data.frames into tables to inspect number of positive, negative and not validated results
validation.DQB1DRB1.m <- table(validation.DQB1DRB1.m$validated)
validation.DQB1DRB1.c <- table(validation.DQB1DRB1.c$validated)
validation.BC.m <- table(validation.BC.m$validated)
validation.BC.c <- table(validation.BC.c$validated)

# Validation plots
# DQB1~DRB1
pie.DQB1DRB1.validation <- data.frame(
  group = c("Positive", "Negative", "Not validated"),
  value = c(sum(validation.DQB1DRB1.m["TRUE"], validation.DQB1DRB1.c["TRUE"])/2,
            sum(validation.DQB1DRB1.m["FALSE"], validation.DQB1DRB1.c["FALSE"])/2,
            sum(validation.DQB1DRB1.m["Not validated"], validation.DQB1DRB1.c["Not validated"])/2)
)

pie.DQB1DRB1.validation <- ggplot(pie.DQB1DRB1.validation, aes(x="", y=value, fill=reorder(group, -value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(fill = "Validation result",
       x = NULL,
       y = NULL,
       title = "HLA-DQB1~DRB1")+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))
pie.DQB1DRB1.validation

# B~C
pie.BC.validation <- data.frame(
  group = c("Positive", "Negative", "Not validated"),
  value = c(sum(validation.BC.m["TRUE"], validation.BC.c["TRUE"])/2,
                       sum(validation.BC.m["FALSE"], validation.BC.c["FALSE"])/2,
                       sum(validation.BC.m["Not validated"], validation.BC.c["Not validated"])/2)
)
pie.BC.validation <- ggplot(pie.BC.validation, aes(x="", y=value, fill=reorder(group, -value)))+
  geom_bar(width = 1, stat = "identity", color="black")+ 
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(fill = "Validation result",
       x = NULL,
       y = NULL,
       title = "HLA-B~C ")+
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5))
pie.BC.validation

# Combined plot with validation results for HLA-B~C and HLA-DQB1~DRB1
grid_arrange_shared_legend(pie.BC.validation, pie.DQB1DRB1.validation, 
             nrow = 2,  
             ncol = 1,
             top = textGrob("Validation of estimated haplotypes", gp=gpar(fontsize=16), vjust = 0.5))
