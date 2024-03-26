library("corrplot")
library("RColorBrewer")

keep <- c("SPONDYLOPATHY_FG", "H7_IRIDOACUTE", "H7_IRIDOCYCLITIS", "K11_REIMB_202", "RX_CROHN_2NDLINE", "RHEUMA_NOS", "M13_POLYARTHROPATHIES", "M13_RHEUMA", "RHEUMA_SEROPOS_OTH")
keep <- unique(melted_df$traits) # melted_df is from hap_regress_plot.R for the heatmap want to use the traits from

# GENETIC COR
gencor <- data.table::fread("/home/ivm/from_satu/general_files/freeze10_combined_gen_cor_sumstats.txt")
gencor$p1 <- gsub(".*munged/","",gencor$p1)
gencor$p1 <- gsub(".sumstats.gz","",gencor$p1)
gencor$p2 <- gsub(".*munged/","",gencor$p2)
gencor$p2 <- gsub(".sumstats.gz","",gencor$p2)

# Filter
gct <- filter(gencor,p1 %in% keep & p2 %in% keep)
  
## Munge and plot gencor
gct <- gct %>% filter(p1!=p2) %>% select(p1,p2,rg,se,p,z) # Drop pairs of traits with themselves
#gct <- gct[!duplicated(t(apply(gct %>% select(p1,p2), 1, sort))),] # get rid of repeat pairs
x2 = reshape2::dcast(gct , p1 ~ p2, value.var = "rg")
mat2 = as.matrix(x2[, 2:ncol(x2)])
rownames(mat2) = x2$p1
mat2[mat2 > 1] = 1
mat2[mat2 < -1] = -1
corr <- mat2
# Use correlation between variables as distance for ordering
dd <- dist(mat2)
hc <- hclust(dd)
mat2[is.na(mat2)] <- 0

# cluster
mat2<-mat2[hc$order, hc$order]
corrplot(mat2,type="lower",tl.col="black",cl.lim = c(-1,1),tl.cex=0.8,order="hclust",is.corr=TRUE, col=colorRampPalette(c("blue","white","red"))(100),method ="color",na.label=" ")

# order same as in heatmap
mat2 <- mat2[match(keep,rownames(mat2)),match(keep,colnames(mat2))]
corrplot(mat2,type="lower",tl.col="black",cl.lim = c(-1,1),tl.cex=0.8,is.corr=TRUE, col=colorRampPalette(c("blue","white","red"))(100),method ="color",na.label=" ")

# PHENOTYPIC COR
p <- data.table::fread("/home/ivm/from_satu/general_files/R10_272_traits_phenotype_corr.gz")
pfilt <- filter(p,Var1 %in% keep & Var2 %in% keep)

## Munge and plot phenotypic cor
x2 = reshape2::dcast(pfilt , Var1 ~ Var2, value.var = "Freq")
mat2 = as.matrix(x2[, 2:ncol(x2)]) # remove trait column
rownames(mat2) = x2$Var1 # add trait list as rownames
mat2[mat2 > 1] = 1
mat2[mat2 < -1] = -1
corr <- mat2
# Use correlation between variables as distance for ordering
dd <- dist(mat2)
hc <- hclust(dd)
mat2[is.na(mat2)] <- 0

# cluster
mat2<-mat2[hc$order, hc$order]
corrplot(mat2,type="lower",tl.col="black",cl.lim = c(-1,1),order="hclust",col=colorRampPalette(c("blue","white","red"))(100),tl.cex=0.8,is.corr=TRUE, method ="color",na.label=" ")

# order same as in heatmap
mat2 <- mat2[match(keep,rownames(mat2)),match(keep,colnames(mat2))]
corrplot(mat2,type="lower",tl.col="black",cl.lim = c(-1,1),col=colorRampPalette(c("blue","white","red"))(100),tl.cex=0.8,is.corr=TRUE, method ="color",na.label=" ")
