## Written by Courtney Smith 10-7-2024
## Goal of script: Permutation analysis to rerun regressions on "random" phenotype assignments

library(dplyr)
library(data.table)
library("RColorBrewer")

# Set parameters
blocks = c(1,2,3)
nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

# PERMUTATION 1: Scramble haplotype group assignments (keeping traits and covariate labels unchanged and correlated with each other as before; note the hap group assignments within a row are still correlated but now the relationship between them and the traits/covariates is scrambled)
for (bnum in blocks){

  print(paste0("Just started: block",bnum," and nsnp",nsnp))

  # Load df with original haplotypes and traits from original regressions
  df <- as.data.frame(data.table::fread(paste0(nsnp_dir,"hb",bnum,"_preregressiondf.txt")))
  trait1_num=which(colnames(df)=="AB1_BACT_INTEST_OTH") # get column id for first trait so can select all trait names
  bin_traits <- colnames(df %>% select(c(trait1_num:ncol(df))))
  keptclusters <- colnames(df %>% select(matches("cluster."),-cluster1a,-cluster1b)) # can use all clusters in df because this already had the to be left out cluster removed
  nclust=length(keptclusters)

  # Scramble IDs
  set.seed(780)
  df[keptclusters] <- lapply(df[keptclusters],sample)

  # Run regression
  results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())

  for (trait in bin_traits) {
    form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=df)
    beta<-summary(regres)$coefficients[2:(nclust+1),1]
    se<-summary(regres)$coefficients[2:(nclust+1),2]
    p<-summary(regres)$coefficients[2:(nclust+1),4]
    tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
    results_dt=rbindlist(list(results_dt,tmp_dt))
  }

  # save regression results
  cluster<-gsub("cluster","",keptclusters)
  results_dt$cluster<-rep(cluster, length(bin_traits))
  data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_regresults_permutation1.txt"), row.names = F, quote = F, sep = "\t")
}

# PERMUTATION 2: Scramble phenotype assignments (keeping hap group assignments and covariate labels unchanged and correlated with each other as before; note the phenotype assignments within a row are still correlated but now the relationship between them and the hap group/covariates is scrambled)
for (bnum in blocks){

  print(paste0("Just started: block",bnum," and nsnp",nsnp))

  # Load df with original haplotypes and traits from original regressions
  df <- as.data.frame(data.table::fread(paste0(nsnp_dir,"hb",bnum,"_preregressiondf.txt")))
  trait1_num=which(colnames(df)=="AB1_BACT_INTEST_OTH") # get column id for first trait so can select all trait names
  bin_traits <- colnames(df %>% select(c(trait1_num:ncol(df))))
  keptclusters <- colnames(df %>% select(matches("cluster."),-cluster1a,-cluster1b)) # can use all clusters in df because this already had the to be left out cluster removed
  nclust=length(keptclusters)

  # Scramble IDs
  set.seed(374)
  df[bin_traits] <- lapply(df[bin_traits],sample)

  # Run regression
  results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())

  for (trait in bin_traits) {
    form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=df)
    beta<-summary(regres)$coefficients[2:(nclust+1),1]
    se<-summary(regres)$coefficients[2:(nclust+1),2]
    p<-summary(regres)$coefficients[2:(nclust+1),4]
    tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
    results_dt=rbindlist(list(results_dt,tmp_dt))
  }

  # save regression results
  cluster<-gsub("cluster","",keptclusters)
  results_dt$cluster<-rep(cluster, length(bin_traits))
  data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_regresults_permutation2.txt"), row.names = F, quote = F, sep = "\t")
}
