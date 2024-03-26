## Courtney Smith - HLA analysis - Cluster regression adjusting for key alleles
## Started September 4th, 2023
## Goal of script: Generate results of cluster regression when adjusting for key alleles, one at a time, and compare to unadjusted cluster regression for all blocks

library(dplyr)
library(data.table)
library(car)
library(ggplot2)

nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

# Adjust all alleles in the respective block
blocks=1:3

for (bnum in blocks){
dat_bin <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_preregressiondf.txt"))
ninfo <- ncol(dat_bin)-269 # only keep columns that have phenotype info, assumes 269 traits
bin_traits <- colnames(dat_bin %>% select(-c(1:ninfo)))
keptclusters <- colnames(dat_bin %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b))
nclust=length(keptclusters)

a <- data.table::fread("/finngen/red/courtrun/R10_amino_allele/R10_alleles98_MAF.txt")
colnames(a) <- gsub("\\*|:","",colnames(a)) # reformat column names to remove special characters
dat_bin <- left_join(dat_bin,a,by=c("FINNGENID"="FID"))

results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())

if (bnum==1){
alleles = colnames(a %>% select(matches("^A")))
}
if (bnum==2){
  alleles = colnames(a %>% select(matches("^B"),matches("^C")))
}
if (bnum==3){
  alleles = colnames(a %>% select(matches("^DR"),matches("^DQ")))
}
justclusters <- keptclusters
keptclusters <- c(keptclusters,alleles)

## First identify those that are independent as defined by all vif < 5
i=0

while (i>(-1)){
  if (i>0){
    keptclusters <- keptclusters[ ! keptclusters %in% c(vifallele)]
  }
  nclust <- length(keptclusters)
  
  # first figure out non colinear ones
  trait="T1D"
  form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
  regres <- glm(form, family="binomial", data=dat_bin)
  vr <- as.data.frame(vif(regres)) %>% mutate(allele=rownames(.)) %>% rename(vif=`vif(regres)`)
  print(nrow(filter(vr,vif>5)))
  if (nrow(filter(vr,vif>5))!=0){
    vifallele=(vr %>% arrange(-vif) %>%head(1))$allele
    i=i+1
  } else{break}
}


## ALL TRAITS
# Run regression
for (trait in bin_traits) {
  form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
  regres <- glm(form, family="binomial", data=dat_bin)
  beta<-summary(regres)$coefficients[2:(nclust+1),1]
  se<-summary(regres)$coefficients[2:(nclust+1),2]
  p<-summary(regres)$coefficients[2:(nclust+1),4]
  tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
  results_dt=rbindlist(list(results_dt,tmp_dt))
  print(trait)
}
print("hi")
# save regression results
results_dt$cluster<-rep(keptclusters, length(bin_traits))

data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_regresults_adj_all_alleles_inblock.txt"), row.names = F, quote = F, sep = "\t")

# Save list of these independent alleles
ja_alleles=data.frame(clustalleles=keptclusters)
data.table::fwrite(ja_alleles, paste0(nsnp_dir,"hb",bnum,"_clustalleles_adj_all_alleles_inblock.txt"), row.names = F, quote = F, sep = "\t")
}

############## method 2 (all alleles in a block together with all haplotypes in the block)
bnum=1
results_dt <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults_adj_all_alleles_inblock.txt"))

# Filter to traits sig in unadjusted cluster regression
cr <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cst <- unique((cr %>% filter(abs(z)>4))$nimi) # keep sig traits

# Label type
clusts <- unique(filter(results_dt,grepl("cluster",cluster))$cluster)
results_dt$type <- ifelse(results_dt$cluster %in% clusts,"cluster","allele")

# Filter to just haplotypes
results_dt_adjallele <- results_dt %>% filter(type=="cluster")
results_dt_adjallele$cluster <- as.integer(gsub("cluster","",results_dt_adjallele$cluster))

# Load in preadjusted reg data
results_dt_justclusters <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt"))

comb <- left_join(results_dt_justclusters %>% mutate(z_pre=beta/se),results_dt_adjallele %>% mutate(z_adj=beta/se),by=c("nimi","cluster"))
ggplot(comb,aes(x=z_pre,y=z_adj))+geom_point(aes(color=as.factor(cluster)))+
  labs(x=paste0("Unadjusted cluster regression for block ",bnum),y=paste0("Cluster regression adjusted for alleles in the block"),color="Cluster")+
  geom_abline(slope=1,intercept=0)+theme_classic()
m <- lm(z_adj ~ z_pre, comb);
format(unname(coef(m)[1]), digits = 2) # y intercept
format(unname(coef(m)[2]), digits = 2) # slope
format(summary(m)$r.squared, digits = 3) # r2

# Number that were sig before that are still sig
nrow(filter(comb,abs(z_pre)>4)) # num sig associations in unadj regression
length(unique((filter(comb,abs(z_pre)>4))$nimi)) # num sig associations in unadj regression
nrow(filter(comb,abs(z_adj)>4)) # num sig associations in adj regression
length(unique((filter(comb,abs(z_adj)>4))$nimi)) # num sig associations in adj regression

head(comb %>% filter(abs(z_pre)>4) %>% arrange(-abs(z_adj)))

bnum=2
results_dt <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults_adj_all_alleles_inblock.txt"))

# Filter to traits sig in unadjusted cluster regression
cr <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cst <- unique((cr %>% filter(abs(z)>4))$nimi) # keep sig traits

# Label type
clusts <- unique(filter(results_dt,grepl("cluster",cluster))$cluster)
results_dt$type <- ifelse(results_dt$cluster %in% clusts,"cluster","allele")

# Filter to just haplotypes
results_dt_adjallele <- results_dt %>% filter(type=="cluster")
results_dt_adjallele$cluster <- as.integer(gsub("cluster","",results_dt_adjallele$cluster))

# Load in preadjusted reg data
results_dt_justclusters <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt"))

comb <- left_join(results_dt_justclusters %>% mutate(z_pre=beta/se),results_dt_adjallele %>% mutate(z_adj=beta/se),by=c("nimi","cluster"))
ggplot(comb,aes(x=z_pre,y=z_adj))+geom_point(aes(color=as.factor(cluster)))+
  labs(x=paste0("Unadjusted cluster regression for block ",bnum),y=paste0("Cluster regression adjusted for alleles in the block"),color="Cluster")+
  geom_abline(slope=1,intercept=0)+theme_classic()
m <- lm(z_adj ~ z_pre, comb);
format(unname(coef(m)[1]), digits = 2) # y intercept
format(unname(coef(m)[2]), digits = 2) # slope
format(summary(m)$r.squared, digits = 3) # r2

# Number that were sig before that are still sig
nrow(filter(comb,abs(z_pre)>4)) # num sig associations in unadj regression
length(unique((filter(comb,abs(z_pre)>4))$nimi)) # num sig associations in unadj regression
nrow(filter(comb,abs(z_adj)>4)) # num sig associations in adj regression
length(unique((filter(comb,abs(z_adj)>4))$nimi)) # num sig associations in adj regression

head(comb %>% filter(abs(z_pre)>4) %>% arrange(-abs(z_adj)))

bnum=3
results_dt <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults_adj_all_alleles_inblock.txt"))

# Filter to traits sig in unadjusted cluster regression
cr <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cst <- unique((cr %>% filter(abs(z)>4))$nimi) # keep sig traits

# Label type
clusts <- unique(filter(results_dt,grepl("cluster",cluster))$cluster)
results_dt$type <- ifelse(results_dt$cluster %in% clusts,"cluster","allele")

# Filter to just haplotypes
results_dt_adjallele <- results_dt %>% filter(type=="cluster")
results_dt_adjallele$cluster <- as.integer(gsub("cluster","",results_dt_adjallele$cluster))

# Load in preadjusted reg data
results_dt_justclusters <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt"))

comb <- left_join(results_dt_justclusters %>% mutate(z_pre=beta/se),results_dt_adjallele %>% mutate(z_adj=beta/se),by=c("nimi","cluster"))
ggplot(comb,aes(x=z_pre,y=z_adj))+geom_point(aes(color=as.factor(cluster)))+
  labs(x=paste0("Unadjusted cluster regression for block ",bnum),y=paste0("Cluster regression adjusted for alleles in the block"),color="Cluster")+
  geom_abline(slope=1,intercept=0)+theme_classic()
m <- lm(z_adj ~ z_pre, comb);
format(unname(coef(m)[1]), digits = 2) # y intercept
format(unname(coef(m)[2]), digits = 2) # slope
format(summary(m)$r.squared, digits = 3) # r2

# Number that were sig before that are still sig
nrow(filter(comb,abs(z_pre)>4)) # num sig associations in unadj regression
length(unique((filter(comb,abs(z_pre)>4))$nimi)) # num sig associations in unadj regression
nrow(filter(comb,abs(z_adj)>4)) # num sig associations in adj regression
length(unique((filter(comb,abs(z_adj)>4))$nimi)) # num sig associations in adj regression

head(comb %>% filter(abs(z_pre)>4) %>% arrange(-abs(z_adj)))

## 

results_dt1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults_adj_all_alleles_inblock.txt")) %>% mutate(block=1)
results_dt2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults_adj_all_alleles_inblock.txt")) %>% mutate(block=2)
results_dt3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults_adj_all_alleles_inblock.txt")) %>% mutate(block=3)
results_dt <- bind_rows(results_dt1,results_dt2,results_dt3)

# Filter to traits sig in unadjusted cluster regression
cr1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cr2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cr3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cr <- bind_rows(cr1,cr2,cr3)
cst <- unique((cr %>% filter(abs(z)>4))$nimi) # keep sig traits

# Label type
clusts <- unique(filter(results_dt,grepl("cluster",cluster))$cluster)
results_dt$type <- ifelse(results_dt$cluster %in% clusts,"cluster","allele")

# Save results of cluster reg adj for all alleles in the block
results_dt_supp <- results_dt %>% mutate(z_adj=beta/se)
data.table::fwrite(results_dt_supp, paste0(nsnp_dir,"regresults_adjallelesinblock_alltraits_supp.txt"), row.names = F, quote = F, sep = "\t")

# Filter to just haplotypes
results_dt_adjallele <- results_dt %>% filter(type=="cluster")
results_dt_adjallele$cluster <- as.integer(gsub("cluster","",results_dt_adjallele$cluster))

# Load in preadjusted reg data
results_dt_justclusters1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults.txt")) %>% mutate(block=1)
results_dt_justclusters2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults.txt")) %>% mutate(block=2)
results_dt_justclusters3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults.txt")) %>% mutate(block=3)
results_dt_justclusters <- bind_rows(results_dt_justclusters1,results_dt_justclusters2,results_dt_justclusters3)

comb <- left_join(results_dt_justclusters %>% mutate(z_pre=beta/se),results_dt_adjallele %>% mutate(z_adj=beta/se),by=c("nimi","cluster","block"))
ggplot(comb,aes(x=z_pre,y=z_adj))+geom_point(aes(color=as.factor(cluster)))+
  labs(x=paste0("Unadjusted cluster regression"),y=paste0("Cluster regression adjusted for alleles in the respective block"),color="Cluster")+
  geom_abline(slope=1,intercept=0)+theme_classic()
m <- lm(z_adj ~ z_pre, comb);
format(unname(coef(m)[1]), digits = 2) # y intercept
format(unname(coef(m)[2]), digits = 2) # slope
format(summary(m)$r.squared, digits = 3) # r2

# Number that were sig before that are still sig
nrow(filter(comb,abs(z_pre)>4)) # num sig associations in unadj regression
length(unique((filter(comb,abs(z_pre)>4))$nimi)) # num sig associations in unadj regression
nrow(filter(comb,abs(z_adj)>4)) # num sig associations in adj regression
length(unique((filter(comb,abs(z_adj)>4))$nimi)) # num sig associations in adj regression
nrow(filter(comb,abs(z_adj)>3)) # num sig associations in adj regression
length(unique((filter(comb,abs(z_adj)>3))$nimi)) # num sig associations in adj regression

head(comb %>% filter(abs(z_pre)>4) %>% arrange(-abs(z_adj)))

# how many were more significant after adjustment across the blocks
(filter(comb,(abs(z_adj)>4|abs(z_pre)>4) & abs(z_adj)>abs(z_pre))) %>% group_by(block) %>% count()
dim(filter(comb,(abs(z_adj)>4|abs(z_pre)>4) & abs(z_adj)>abs(z_pre)))


######## chi squared analysis on method 2 (all alleles in a block together with all haplotypes in the block)
bnum=1
results_dt <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults_adj_all_alleles_inblock.txt"))

# Filter to traits sig in unadjusted cluster regression
cr <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cst <- unique((cr %>% filter(abs(z)>4))$nimi) # keep sig traits

# Label type
clusts <- unique(filter(results_dt,grepl("cluster",cluster))$cluster)
results_dt$type <- ifelse(results_dt$cluster %in% clusts,"cluster","allele")

# Calculate the expected quantiles under the null hypothesis (uniform distribution)
results_dt <- results_dt %>% filter(nimi %in% cst) %>% mutate(o=-log10(p)) %>% arrange(o) %>%
  group_by(nimi,type) %>% mutate(e=-log10(rank(p)/n()),e2=-log10(ppoints(p)))

library(ggplot2)

ggplot(results_dt %>% filter(nimi=="AUTOIMMUNE"),aes(x=e,y=o,color=type)) + facet_wrap(~nimi)+
  geom_point()+geom_abline(intercept=0,slope=1)+
  labs(x="Expected -log10(P)",y="Observed -log10(P)")+
  theme_classic()

ggplot(results_dt %>% filter(nimi=="K11_OTHDIG"),aes(x=e,y=o,color=type)) + facet_wrap(~nimi)+
  geom_point()+geom_abline(intercept=0,slope=1)+
  labs(x="Expected -log10(P)",y="Observed -log10(P)")+
  theme_classic()

ggplot(results_dt,aes(x=e,y=o,color=type)) + facet_wrap(~nimi)+
  geom_point()+geom_abline(intercept=0,slope=1)+
  labs(x="Expected -log10(P)",y="Observed -log10(P)",color="Predictor Type")+
  theme_classic()
