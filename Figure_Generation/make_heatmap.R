# Courtney Smith - HLA Project - HLA Heatmap Plot
# Script started 8-3-2022
# Goal of script: Plot heatmap of zscores for traits vs variants ascertained in them; with and without x axis

library("ggplot2")
library("ggdendro")
library("reshape2")
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggbiplot)
library(forcats)
library("pheatmap")

full <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R8_sumstats/matrix/filt185hits_wide.txt") # commandArgs(TRUE)[1]
t <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/hla/uncorr185traits/uncorr_traits185.txt",header=F) # commandArgs(TRUE)[2]
full <- full %>% select(ID,matches(t$V1)) # subset columns to just rsids and traits in t
extra <- setdiff(gsub("beta_","",colnames(full %>% select(matches("^beta_")))),t$V1) # identify and then remove extra traits that got included that shouldn't have
full <- full %>% select(-matches(extra),matches("M13_ARTHROSIS_KNEE_PRIM_KNEESURG"))

betas_full <- full %>% select(matches("^beta_"))
colnames(betas_full) <- gsub("beta_","",colnames(betas_full))
se_full <- full %>% select(matches("sebeta_"))
z <- betas_full/se_full
z <- t(z)
my_function <- function(x) {return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}
z <- t(apply(z,FUN=my_function,MARGIN=1)) # apply transform to each row (1=row, 2=column), row normalize
means <- rowMeans(as.matrix(z)) # calculate the mean of each row (aka each Trait) of z
SDs <- rowSds(as.matrix(z)) # calculate the SD of each row (aka each Trait) of z
z <- sweep(z, 1, means) # subtract each row by its respective mean
z <- sweep(z, 1, SDs, FUN = '/') # and divide each row by its respective SD
colnames(z) <- full$ID
z <- na.omit(z)
z <- t(z)
z[rowMedians(as.matrix(z))<0,] <- z[rowMedians(as.matrix(z))<0,]*-1 # flip sign so rsids have median zscore across Traits thats positive
pvalues <- full %>% select(ID,matches("pval_"))
z_full <- as.data.frame(z)
z_full$ID <- rownames(z_full)
#z.long <- melt(z, id = c("Trait"))

# Add coloring to Trait labels
an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/hla/uncorr185traits/traits_new_class.txt") %>% select(trait,class) # commandArgs(TRUE)[3]
an <- as.data.frame(an)
an <- right_join(an,data.frame(trait=colnames(z))) %>% arrange(trait)
rownames(an) <- an$trait
an <- an %>% select(class)

# Add labels for genes
ang <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/hla/uncorr185traits/genes_boundaries_rsids_an.txt")
ang <- filter(ang,ID %in% rownames(z)) %>% arrange(pos)
ang <- left_join(data.frame(ID=rownames(z)),ang,by=c("ID")) %>% select(ID,nearest_genes)
ang <- ang %>% mutate(label = paste0(ID,"-",nearest_genes))
filter(data.frame(V1=ang$ID,V2=rownames(z)),V1!=V2) # make sure they are lines up properly (should be empty)
rownames(z) <- ang$label

png("/oak/stanford/groups/pritch/users/strausz/finngen_R8_sumstats/figures/pheat_185fullhits_alltraits.png",width=4000,height=4000) # commandArgs(TRUE)[4])
pheatmap(t(z), annotation_row=an, show_colnames=T, cluster_cols=T, cluster_rows=T)
dev.off()

ang2 <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/hla/uncorr185traits/genes_boundaries_rsids_an.txt")
zsnp <- left_join(z_full,distinct(ang2 %>% select(ID,pos)),by=c("ID")) %>% arrange(pos)
zsnp <- zsnp %>% left_join(ang,by=c("ID"))
rownames(zsnp) <- zsnp$label
zsnp <- zsnp %>% select(-ID,-nearest_genes,-label,-pos)

png("/oak/stanford/groups/pritch/users/strausz/finngen_R8_sumstats/figures/pheat_185fullhits_alltraits_pos.png",width=4000,height=4000) # commandArgs(TRUE)[4])
pheatmap(t(zsnp), annotation_row=an, show_colnames=T, cluster_rows=T, cluster_cols=F)
dev.off()

# SAVE z matrix
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R8_sumstats/matrix/filt185fullhits_wide_stdz.txt"
write.table(z_full, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)
