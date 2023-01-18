# Courtney Smith - HLA Analysis - LDSC Genetic Correlation
# Goal of script: Plot LDSC genetic correlation matrix

# ml R/4.0
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(corrplot)

combined_genetic_correlations=commandArgs(TRUE)[1] # "/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/combined_gen_cor_sumstats.txt"
gencor <- data.table::fread(combined_genetic_correlations,fill=TRUE)

# get rid of pairs where one or more files did not exist
gencor <- na.omit(gencor)

# Reformat trait columns to just be metabolite names
gencor$p1 <- gsub(".*munged/","",gencor$p1)
gencor$p1 <- gsub(".sumstats.gz","",gencor$p1)
gencor$p2 <- gsub(".*munged/","",gencor$p2)
gencor$p2 <- gsub(".sumstats.gz","",gencor$p2)

# Filter to HLA traits
t <- data.table::fread(commandArgs(TRUE)[2],header=F) # "/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/independent_traits_list.tsv"
gencor <- filter(gencor,p1 %in% t$V1 & p2 %in% t$V1)

# Switch this into wide matrix form
x2 = dcast(gencor , p1 ~ p2, value.var = "rg")
mat2 = as.matrix(x2[, 2:ncol(x2)])
rownames(mat2) = x2$p1
mat2[mat2 > 1] = 1
mat2[mat2 < -1] = -1
corr <- mat2

# Use correlation between variables as distance for ordering
dd <- dist(mat2)
hc <- hclust(dd)
mat2<-mat2[hc$order, hc$order]
mat2[upper.tri(mat2,diag=TRUE)]<- NA # half to NAs

cm <- as.data.frame(mat2) %>% mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>%
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE)

corr <- as.data.frame(as.table(corr))
corr <- na.omit(corr) #remove the NA values
#sort by highest correlation
corr <- corr[order(-abs(corr$Freq)),]
mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
mtx_corr[is.na(mtx_corr)] <- 0

# Add coloring to trait labels
cat <- data.table::fread(commandArgs(TRUE)[3]) %>% select(trait,Category) # "/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/traits_annotation.tsv"

#plot correlations visually
pdf(commandArgs(TRUE)[4],height=30,width=28) # "/oak/stanford/groups/pritch/users/courtrun/projects/hla/figures/LDSC.gencorrplot.corrplot.pdf"
corrplot(mtx_corr,type="lower",tl.col="black",col.lim = c(-1,1),order="hclust",is.corr=FALSE, method = c("square"),na.label=" ")
dev.off()

### Make separate plots per category
for (categ in unique(cat$Category)){
keep_trait <- filter(cat,Category %in% categ)$trait
gencorfilt <- filter(gencor,p1 %in% keep_trait & p2 %in% keep_trait)

# Switch this into wide matrix form
x2 = dcast(gencorfilt, p1 ~ p2, value.var = "rg")
mat2 = as.matrix(x2[, 2:ncol(x2)])
rownames(mat2) = x2$p1
mat2[mat2 > 1] = 1
mat2[mat2 < -1] = -1
corr <- mat2

# Use correlation between variables as distance for ordering
dd <- dist(mat2)
hc <- hclust(dd)
mat2<-mat2[hc$order, hc$order]
mat2[upper.tri(mat2,diag=TRUE)]<- NA # half to NAs

cm <- as.data.frame(mat2) %>% mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>%
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE)

corr <- as.data.frame(as.table(corr))
corr <- na.omit(corr) #remove the NA values
#sort by highest correlation
corr <- corr[order(-abs(corr$Freq)),]
mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
mtx_corr[is.na(mtx_corr)] <- 0

#plot correlations visually
pdf(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/figures/LDSC.gencorrplot.corrplot_",categ,".pdf"),height=20,width=18)
corrplot(mtx_corr,type="lower",tl.col="black",col.lim = c(-1,1),order="hclust",is.corr=FALSE, method = c("square"),na.label=" ")
dev.off()
}
