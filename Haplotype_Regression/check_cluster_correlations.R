# Check correlations between clusters across blocks
library(dplyr)

nsnp=1000

hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

b1 <- data.table::fread("hb1_clusterdoses.txt")
b2 <- data.table::fread("hb2_clusterdoses.txt")
b3 <- data.table::fread("hb3_clusterdoses.txt")

temp <- data.frame()

for (ncol1 in 2:ncol(b1)) {
  for (ncol2 in 2:ncol(b2)) {
    result=cor.test(b1[[ncol1]],b2[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=ncol2-1,corr=corr,corp=corp,blocka=1,blockb=2))
  }
}

for (ncol1 in 2:ncol(b1)) {
  for (ncol2 in 2:ncol(b3)) {
    result=cor.test(b1[[ncol1]],b3[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=ncol2-1,corr=corr,corp=corp,blocka=1,blockb=3))
  }
}

for (ncol1 in 2:ncol(b2)) {
  for (ncol2 in 2:ncol(b3)) {
    result=cor.test(b2[[ncol1]],b3[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=ncol2-1,corr=corr,corp=corp,blocka=2,blockb=3))
  }
}

cor_results <- bind_rows(temp)
data.table::fwrite(cor_results, paste0(nsnp_dir,"correlation_acrossblocks.txt"), row.names = F, quote = F, sep = "\t")

cor_results <- data.table::fread(paste0(nsnp_dir,"correlation_acrossblocks.txt"))
head(cor_results %>% arrange(-abs(corr)),35)
summary(abs(cor_results$corr))

#### Correlations within a block

temp <- data.frame()

for (ncol1 in 2:ncol(b1)) {
  for (ncol2 in 2:ncol(b1)) {
    result=cor.test(b1[[ncol1]],b1[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=ncol2-1,corr=corr,corp=corp,blocka=1,blockb=1))
  }
}

for (ncol1 in 2:ncol(b2)) {
  for (ncol2 in 2:ncol(b2)) {
    result=cor.test(b2[[ncol1]],b2[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=ncol2-1,corr=corr,corp=corp,blocka=2,blockb=2))
  }
}


for (ncol1 in 2:ncol(b3)) {
  for (ncol2 in 2:ncol(b3)) {
    result=cor.test(b3[[ncol1]],b3[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=ncol2-1,corr=corr,corp=corp,blocka=3,blockb=3))
  }
}

cor_results <- bind_rows(temp)
dim(filter(cor_results,Cluster_a!=Cluster_b & abs(corr)>0.5)) # make sure to count in half since duplicate with reverse pair
data.table::fwrite(cor_results, paste0(nsnp_dir,"correlation_withinblocks.txt"), row.names = F, quote = F, sep = "\t")

cor_results <- data.table::fread(paste0(nsnp_dir,"correlation_withinblocks.txt")) %>% filter(Cluster_a!=Cluster_b) 
nline=76
head(cor_results %>% arrange(-abs(corr)),nline)[rep(c(TRUE,FALSE),nline/2),] # don't print the same thing in reverse
summary(abs(((cor_results %>% arrange(-abs(corr)))[rep(c(TRUE,FALSE),nrow(cor_results)/2),])$corr))

#### Between clusters and alleles
a <- data.table::fread("/finngen/red/courtrun/R10_amino_allele/R10_alleles98_MAF.txt")

temp <- data.frame()

for (ncol1 in 2:ncol(b1)) {
  for (ncol2 in 2:ncol(a)) {
    result=cor.test(b1[[ncol1]],a[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    allelename=colnames(a)[ncol2]
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=allelename,corr=corr,corp=corp,blocka=1,blockb="allele"))
  }
}

for (ncol1 in 2:ncol(b2)) {
  for (ncol2 in 2:ncol(a)) {
    result=cor.test(b2[[ncol1]],a[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    allelename=colnames(a)[ncol2]
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=allelename,corr=corr,corp=corp,blocka=2,blockb="allele"))
  }
}

for (ncol1 in 2:ncol(b3)) {
  for (ncol2 in 2:ncol(a)) {
    result=cor.test(b3[[ncol1]],a[[ncol2]],method="spearman",exact=FALSE)
    corr=result$estimate
    corp=result$p.value
    allelename=colnames(a)[ncol2]
    temp <- rbind(temp,data.frame(Cluster_a=ncol1-1,Cluster_b=allelename,corr=corr,corp=corp,blocka=3,blockb="allele"))
  }
}

cor_results <- bind_rows(temp)
data.table::fwrite(cor_results, paste0(nsnp_dir,"correlation_withalleles.txt"), row.names = F, quote = F, sep = "\t")

cor_results <- data.table::fread(paste0(nsnp_dir,"correlation_withalleles.txt"))
summary(abs(cor_results$corr))
