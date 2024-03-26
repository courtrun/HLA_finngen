library(dplyr)
library("corrplot")
library("RColorBrewer")

nsnp=1000

hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

gct <- data.table::fread(paste0(nsnp_dir,"correlation_withalleles.txt")) %>%
  mutate(Clusterblock=paste0("b",blocka,"c",Cluster_a)) %>%
  rename(p1=Clusterblock,p2=Cluster_b)

## Munge and plot gencor
x2 = reshape2::dcast(gct , p1 ~ p2, value.var = "corr")
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
#mat2<-mat2[hc$order, hc$order]
corrplot(mat2,is.corr=FALSE,tl.cex=0.6,cl.lim=c(-1,1),tl.col="black",col=colorRampPalette(c("blue","white","red"))(100),method ="color",na.label=" ")

gct <- gct %>% mutate(r2=corr^2)
filter(gct,r2>0.8)
