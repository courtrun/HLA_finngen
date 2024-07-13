library(dplyr)
library(ggplot2)

library(dplyr)
library(data.table)
library(car)
library(ggplot2)

nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

dat_bin <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_preregressiondf.txt"))
ninfo <- ncol(dat_bin)-269 # only keep columns that have phenotype info, assumes 269 traits
bin_traits <- colnames(dat_bin %>% select(-c(1:ninfo)))
keptclusters <- colnames(dat_bin %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b))
nclust=length(keptclusters)

bnum=1

### Traits to consider
#lupus/RA and such with each other, with IBD
#SPONDYLOPATHY_FG
#L12_LUPUS
#K11_IBD_STRICT
#K11_IBD_STRICT_PSC
#RHEUMA_NOS#

#AB1_EBV and G6_MS

#RX_CROHN_1STLINE and IBD / RHEUMA / infections

#infection: AB1_ERYSIPELAS, AB1_OTHER_SEPSIS, AB1_SEXUAL_TRANSMISSION
#cancer: C3_SQUOMOUS_CELL_CARCINOMA_SKIN_EXALLC, C3_BASAL_CELL_CARCINOMA_EXALLC, 

#T1D and eye things (H7_GLAUCOMA, H7_IRIDOACUTE, H7_IRIDOCHRONIC, H7_CATARACTOTHER)
#cancer and diabetes/infection
#lupus/ RA and vasculitis (I9_ATHSCLE,FG_DOAAC)

trait1="K11_IBD_STRICT"
trait2="RX_CROHN_2NDLINE"

form <- as.formula(paste(paste(trait1,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+",trait2)))
regres <- glm(form, family="binomial", data=dat_bin)
beta<-summary(regres)$coefficients[2:(nclust+1),1]
se<-summary(regres)$coefficients[2:(nclust+1),2]
p<-summary(regres)$coefficients[2:(nclust+1),4]
tmp_dt=data.table("nimi"=trait1,"beta"=beta,"se"=se,"p"=p,trait2=trait2)

form <- as.formula(paste(paste(trait1,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
regres <- glm(form, family="binomial", data=dat_bin)
beta<-summary(regres)$coefficients[2:(nclust+1),1]
se<-summary(regres)$coefficients[2:(nclust+1),2]
p<-summary(regres)$coefficients[2:(nclust+1),4]
tmp_dt_unadj=data.table("nimi"=trait1,"beta"=beta,"se"=se,"p"=p)

df <- bind_cols(tmp_dt,tmp_dt_unadj,by=c("nimi"))
df[,c(2,4,7,9)] %>% filter(`p...4`<0.05|`p...9`<0.05)

ggplot(df,aes(x=-log10(`p...9`),y=-log10(`p...4`)))+geom_point()+
  geom_abline(slope=1,intercept=0)+
  labs(x=paste0("-Log10 P-values for Cluster Associations for Block 1\nfor Unadjusted ",trait1),
       y=paste0("-Log10 P-values for Cluster Associations for Block 1\nfor ",trait1, " Adjusted for ",trait2)) + theme_classic()
