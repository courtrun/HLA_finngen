## Courtney Smith - HLA analysis - Allele regression
## Started September 6th, 2023
## Goal of script: Generate results of allele regression (two ways: one with all vif < 5 independent alleles at once and one, one with one allele at a time)

library(dplyr)
library(data.table)
library(car)
library(ggplot2)

nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

######### RUN Regression just with each allele with each trait directly BUT only one allele at a time in each regression on
## For all traits in finngen!
nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

a <- data.table::fread("/finngen/red/courtrun/R10_amino_allele/R10_alleles98_MAF.txt")
colnames(a) <- gsub("\\*|:","",colnames(a)) # reformat column names to remove special characters

# Get traits
cov<-fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_PHENO_V1.txt.gz")
traits <- fread("/finngen/green/courtrun/finngen_R10_endpoint_core_noncore_1.0.txt") %>%
  filter(CORE_ENDPOINTS=="yes") %>% select(NAME) %>% rename(trait=NAME)
traits<-traits[!duplicated(traits$trait),]
traits<-traits[,c(trait)]
rest<-c("FINNGENID","AGE_AT_DEATH_OR_END_OF_FOLLOWUP","SEX_IMPUTED","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
vars<-append(rest, traits)
dat<-cov[,vars, with=F]
dat[is.na(dat)]<-0
bin_traits <- traits

# Df with alleles and covariates
dat<-merge(a,dat, by.y = "FINNGENID", by.x = "FID")

results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric(),allele=character())

keptclusters <- colnames(a %>% select(-FID)) # list of alleles to run regression for
nclust <- 1

for (trait in bin_traits) {
  for (allele in keptclusters){
    form <- as.formula(paste(paste(trait,"~ "),allele,paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=dat)
    beta<-summary(regres)$coefficients[2:(nclust+1),1]
    se<-summary(regres)$coefficients[2:(nclust+1),2]
    p<-summary(regres)$coefficients[2:(nclust+1),4]
    tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p,allele=allele)
    results_dt=rbindlist(list(results_dt,tmp_dt))
    print(allele)
  }
  print(trait)
}

data.table::fwrite(results_dt, paste0(nsnp_dir,"alleles_indiv_regresults_alltraits.txt"), row.names = F, quote = F, sep = "\t")

### Run on all independent alleles together on all traits - broken out by block
a <- data.table::fread("/finngen/red/courtrun/R10_amino_allele/R10_alleles98_MAF.txt")
colnames(a) <- gsub("\\*|:","",colnames(a)) # reformat column names to remove special characters

# Get traits
cov<-fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_PHENO_V1.txt.gz")
traits <- fread("/finngen/green/courtrun/finngen_R10_endpoint_core_noncore_1.0.txt") %>%
  filter(CORE_ENDPOINTS=="yes") %>% select(NAME) %>% rename(trait=NAME)
traits<-traits[!duplicated(traits$trait),]
traits<-traits[,c(trait)]
rest<-c("FINNGENID","AGE_AT_DEATH_OR_END_OF_FOLLOWUP","SEX_IMPUTED","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
vars<-append(rest, traits)
dat<-cov[,vars, with=F]
dat[is.na(dat)]<-0
bin_traits <- traits

# Df with alleles and covariates
dat<-merge(a,dat, by.y = "FINNGENID", by.x = "FID")
## First identify alleles that are independent as defined by all vif < 5
keptclusters <- colnames(a %>% select(-FID)) # list of alleles to run regression for
i=0

allele_blockA_list <- keptclusters[c(1:10)] # all A alleles (block 1)
allele_blockB_list <- keptclusters[c(11:42)] # all B and C alleles (block 2)
allele_blockC_list <- keptclusters[c(43:98)] # all D* alleles (block 3)

block_lists <- list(c(allele_blockA_list),c(allele_blockB_list),c(allele_blockC_list))

for (block_list in block_lists){
keptclusters <- unlist(block_list)
block_letter <- substring(keptclusters[1],1,1)
while (i>(-1)){
  if (i>0){
    keptclusters <- keptclusters[ ! keptclusters %in% c(vifallele)]
  }
  nclust <- length(keptclusters)
  
  # first figure out non colinear ones
  trait="T1D"
  form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
  regres <- glm(form, family="binomial", data=dat)
  vr <- as.data.frame(vif(regres)) %>% mutate(allele=rownames(.)) %>% rename(vif=`vif(regres)`)
  print(nrow(filter(vr,vif>5)))
  if (nrow(filter(vr,vif>5))!=0){
    vifallele=(vr %>% arrange(-vif) %>%head(1))$allele
    i=i+1
  } else{break}
}

### Run regression on all traits for these remaining independent alleles

results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())

for (trait in bin_traits) {
  form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
  regres <- glm(form, family="binomial", data=dat)
  beta<-summary(regres)$coefficients[2:(nclust+1),1]
  se<-summary(regres)$coefficients[2:(nclust+1),2]
  p<-summary(regres)$coefficients[2:(nclust+1),4]
  tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
  results_dt=rbindlist(list(results_dt,tmp_dt))
}

# save regression results
results_dt$allele<-rep(keptclusters, length(bin_traits))
data.table::fwrite(results_dt, paste0(nsnp_dir,"alleles_vifindep_block",block_letter,"_regresults_alltraits.txt"), row.names = F, quote = F, sep = "\t")
}
