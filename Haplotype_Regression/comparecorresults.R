library(dplyr)
library(ggplot2)
library(data.table)

nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

# Get traits sig in at least one haplotype block regression
f1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regdatatoplot.txt"),header=T)$V1
f2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regdatatoplot.txt"),header=T)$V1
f3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regdatatoplot.txt"),header=T)$V1
keep <- unique(c(f1,f2,f3))

## Block 1
bnum=1
f <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regdatafull.txt"),header=T) %>% filter(V1 %in% keep)
ft <- t(f %>% select(-V1))
colnames(ft) <- f$V1
df_full <- data.frame()
for (t1 in 1:ncol(ft)){
for (t2 in 1:ncol(ft)){
  s1 <- ft[,t1]
  s2 <- ft[,t2]
  p <- cor.test(s1,s2,method="pearson")
  sp <- cor.test(s1,s2,method="spearman")
  df <- data.frame(pearson_cor=p$estimate,
                   selow=p$estimate-(p$estimate-p$conf.int[1])/1.96,
                   sehigh=p$estimate+(-p$estimate+p$conf.int[2])/1.96,
                   pearson_p=p$p.value,
                   spearman_cor=sp$estimate,
                   spearman_p=sp$p.value) %>% mutate(trait1=colnames(ft)[t1],trait2=colnames(ft)[t2])
  df_full=rbindlist(list(df_full,df))
}
}
df_full1 <- df_full %>% filter(trait1!=trait2) %>% 
  filter(trait1 %in% keep & trait2 %in% keep) %>%
  mutate(traitpair=paste0(trait1,"-",trait2))

# Block 2
bnum=2
f <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regdatafull.txt"),header=T)  %>% filter(V1 %in% keep)
ft <- t(f %>% select(-V1))
colnames(ft) <- f$V1
df_full <- data.frame()
for (t1 in 1:ncol(ft)){
  for (t2 in 1:ncol(ft)){
    s1 <- ft[,t1]
    s2 <- ft[,t2]
    p <- cor.test(s1,s2,method="pearson")
    sp <- cor.test(s1,s2,method="spearman")
    df <- data.frame(pearson_cor=p$estimate,
                     selow=p$estimate-(p$estimate-p$conf.int[1])/1.96,
                     sehigh=p$estimate+(-p$estimate+p$conf.int[2])/1.96,
                     pearson_p=p$p.value,spearman_cor=sp$estimate,spearman_p=sp$p.value) %>% mutate(trait1=colnames(ft)[t1],trait2=colnames(ft)[t2])
    df_full=rbindlist(list(df_full,df))
  }
}
df_full2 <- df_full %>% filter(trait1!=trait2) %>% filter(trait1 %in% keep & trait2 %in% keep) %>%
  mutate(traitpair=paste0(trait1,"-",trait2))

## Block 3
bnum=3
f <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regdatafull.txt"),header=T)  %>% filter(V1 %in% keep)
ft <- t(f %>% select(-V1))
colnames(ft) <- f$V1
df_full <- data.frame()
for (t1 in 1:ncol(ft)){
  for (t2 in 1:ncol(ft)){
    s1 <- ft[,t1]
    s2 <- ft[,t2]
    p <- cor.test(s1,s2,method="pearson")
    sp <- cor.test(s1,s2,method="spearman")
    df <- data.frame(pearson_cor=p$estimate,
                     selow=p$estimate-(p$estimate-p$conf.int[1])/1.96,
                     sehigh=p$estimate+(-p$estimate+p$conf.int[2])/1.96,
                     pearson_p=p$p.value,spearman_cor=sp$estimate,spearman_p=sp$p.value) %>% mutate(trait1=colnames(ft)[t1],trait2=colnames(ft)[t2])
    df_full=rbindlist(list(df_full,df))
  }
}
df_full3 <- df_full %>% filter(trait1!=trait2) %>% filter(trait1 %in% keep & trait2 %in% keep) %>%
  mutate(traitpair=paste0(trait1,"-",trait2))

# GENETIC COR
gencor <- data.table::fread("/home/ivm/from_satu/general_files/freeze10_combined_gen_cor_sumstats.txt")
gencor$p1 <- gsub(".*munged/","",gencor$p1)
gencor$p1 <- gsub(".sumstats.gz","",gencor$p1)
gencor$p2 <- gsub(".*munged/","",gencor$p2)
gencor$p2 <- gsub(".sumstats.gz","",gencor$p2)
gct <- filter(gencor,p1 %in% keep & p2 %in% keep)
gct <- gct %>% filter(p1!=p2) %>% select(p1,p2,rg,se,p,z) %>% 
  mutate(traitpair=paste0(p1,"-",p2),selow=rg-se,sehigh=rg+se)

## PHENOTYPIC COR
p <- data.table::fread("/home/ivm/from_satu/general_files/R10_272_traits_phenotype_corr.gz")
se <- reshape2::melt(data.table::fread("/home/ivm/from_satu/general_files/R10_pheno_corr/R10_pheno_corr_SE.txt"))
pval <- reshape2::melt(data.table::fread("/home/ivm/from_satu/general_files/R10_pheno_corr/R10_pheno_corr_pval.txt"))
pfilt <- filter(p,Var1 %in% keep & Var2 %in% keep) %>% mutate(traitpair=paste0(Var2,"-",Var1))
pfilt <- left_join(left_join(pfilt,se,by=c("Var1"="V1","Var2"="variable")),pval,by=c("Var1"="V1","Var2"="variable")) %>%
  rename(se=value.x,p=value.y,phenocor=Freq) %>% mutate(selow=phenocor-se,sehigh=phenocor+se)

## Combine
cw <- left_join(left_join(left_join(left_join(df_full1,gct,by=c("traitpair")),pfilt,by=c("traitpair")),df_full2,by=c("traitpair")),df_full3,by=c("traitpair"))
cl <- bind_rows(df_full1 %>% select(pearson_cor,traitpair,selow,sehigh) %>% rename(cor=pearson_cor) %>% mutate(measure="hapregb1"),
                df_full2 %>% select(pearson_cor,traitpair,selow,sehigh) %>% rename(cor=pearson_cor) %>% mutate(measure="hapregb2"),
                df_full3 %>% select(pearson_cor,traitpair,selow,sehigh) %>% rename(cor=pearson_cor) %>% mutate(measure="hapregb3"),
                gct %>% select(rg,traitpair,selow,sehigh) %>% rename(cor=rg) %>% mutate(measure="gc"),
                pfilt %>% select(phenocor,traitpair,selow,sehigh) %>% rename(cor=phenocor) %>% mutate(measure="pc"))
cl$cor <- as.numeric(ifelse(cl$cor>1,1,ifelse(cl$cor<(-1),-1,cl$cor)))

ggplot(cl,aes(x=measure,y=reorder(traitpair,cor),fill=cor))+
  geom_tile()+labs(x="Correlation Measures",y="Trait Pairs",fill="Correlation")+
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.text=element_text(size=20),
        axis.title = element_text(size=24),
        legend.text=element_text(size=20),
        legend.title=element_text(size=24))+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits = c(-1, 1))+
  scale_x_discrete(labels=c("Genetic\nCorrelation","Hap Reg\nBlock 1","Hap Reg\nBlock 2","Hap Reg\nBlock 3","Phenotypic\nCorrelation"))

sharedtraits <- allsig # from hap_regress_plot_analysis.R
t <- cl %>% group_by(traitpair) %>% #filter(traitpair %in% sharedtraits) %>%
  filter(any(cor > 0.8&measure=="hapregb1"))%>% 
  filter(any(cor > 0.8&measure=="hapregb2"))%>%
  filter(any(cor > 0.8&measure=="hapregb3")) %>% arrange(traitpair)

#################
melted_df <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_heatmapplotted.txt"))

cl %>% filter(abs(cor) > 0.9) %>% group_by(traitpair) %>% count() %>%
  filter(n>3)
cl %>% filter(traitpair=="H7_IRIDOCYCLITIS-H7_IRIDOACUTE")
filter(melted_df,grepl("IRIDO",traits))

cl %>% filter(traitpair=="E4_GRAVES_STRICT-M13_RHEUMA")

nrow(filter(cl,cor>0.3&measure=="gc"))
nrow(filter(cl,cor<0.3&measure=="gc"))
nrow(filter(cl,cor<(-0.3)&measure=="gc"))

length(unique((cl %>% group_by(traitpair) %>%
  filter(any(cor>0.3&measure=="gc")&
           any(cor>0.3&measure=="hapregb1")&
         any(cor>0.3&measure=="hapregb2")&
           any(cor>0.3&measure=="hapregb3")))$traitpair))

length(unique((cl %>% group_by(traitpair) %>%
                 filter(any(cor<(-0.3)&measure=="gc")&
                          (any(cor>0.3&measure=="hapregb1")|
                          any(cor>0.3&measure=="hapregb2")|
                          any(cor>0.3&measure=="hapregb3"))))$traitpair))

length(unique((cl %>% group_by(traitpair) %>%
                 filter(any(cor>0.3&measure=="gc")&(
                          any(cor<(-0.3)&measure=="hapregb1")|
                          any(cor<(-0.3)&measure=="hapregb2")|
                          any(cor<(-0.3)&measure=="hapregb3"))))$traitpair))

length(unique((cl %>% group_by(traitpair) %>%
                 filter(any(cor>0.3&measure=="gc")&(
                   any(cor<(-0.3)&measure=="hapregb1")&
                     any(cor<(-0.3)&measure=="hapregb2")&
                     any(cor<(-0.3)&measure=="hapregb3"))))$traitpair))

length(unique((cl %>% group_by(traitpair) %>%
                 filter(any(cor<0.3&measure=="gc")&
                          any(cor>0.3&measure=="hapregb1")&
                          any(cor>0.3&measure=="hapregb2")&
                          any(cor>0.3&measure=="hapregb3")))$traitpair))

######################

filttraits <- unique((cl %>% filter(cor<(-0.6)))$traitpair)
ggplot(cl %>% filter(traitpair %in% filttraits),aes(x=measure,y=reorder(traitpair,cor),fill=cor))+
  geom_tile()+labs(x="Correlation Measures",y="Trait Pairs")+
  theme_classic()+
  #theme(axis.text.y=element_blank())+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits = c(-1, 1))+
  scale_x_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))

## Pick out interesting trait pairs
unique((filter(cw,pearson_cor.x<(-0.5)&rg>0.4)) %>% select(trait1,trait2))
filter(cw,traitpair %in% c("G6_MS-AUTOIMMUNE_NONTHYROID",
"G6_MS-RX_CROHN_2NDLINE",
"G6_OTHDEMYEL-RHEUMA_SEROPOS_OTH"))

unique((filter(cw,pearson_cor>(0.3)&rg<(-0.3))) %>% select(trait1,trait2,rg,pearson_cor))
filter(cw,traitpair %in% c("AB1_SEXUAL_TRANSMISSION_NOS-L12_LUPUS",
       "AB1_SEXUAL_TRANSMISSION_NOS-N14_GLOMEINOTH"))

unique((filter(cw,pearson_cor<(-0.5)&trait1=="E4_GRAVES_STRICT")) %>% select(trait1,trait2,rg,pearson_cor))
filter(cw,traitpair %in% c("E4_GRAVES_STRICT-M13_RHEUMA",
                           "E4_GRAVES_STRICT-M13_POLYARTHROPATHIES",
                           "E4_GRAVES_STRICT-RHEUMA_NOS"))

unique((filter(cw,phenocor>(0.3)&pearson_cor<(0.2))) %>% select(trait1,trait2,rg,pearson_cor))
filter(cw,traitpair %in% c("AUTOIMMUNE_NONTHYROID-K11_ENERCOLNONINF",
"AUTOIMMUNE_NONTHYROID-K11_IBD_STRICT"))

unique((filter(cw,rg>(0.7)&pearson_cor<(-0.2))) %>% select(trait1,trait2,rg,pearson_cor))
filter(cw,traitpair=="K11_ENERCOLNONINF-M13_SACROILIITIS")

unique((filter(cw,rg>(0.5)&pearson_cor<(-0.5))) %>% select(trait1,trait2,rg,pearson_cor))
filter(cw,traitpair=="AUTOIMMUNE_HYPERTHYROIDISM-RHEUMA_SEROPOS_OTH")

tps <- c("G6_MS-AUTOIMMUNE_NONTHYROID","AB1_SEXUAL_TRANSMISSION_NOS-L12_LUPUS",
         "E4_GRAVES_STRICT-M13_RHEUMA","AUTOIMMUNE_NONTHYROID-K11_ENERCOLNONINF",
         "K11_ENERCOLNONINF-M13_SACROILIITIS",
         "AUTOIMMUNE_HYPERTHYROIDISM-RHEUMA_SEROPOS_OTH",
         "G6_MS-RX_CROHN_2NDLINE",
         "G6_OTHDEMYEL-RHEUMA_SEROPOS_OTH",
         "AB1_SEXUAL_TRANSMISSION_NOS-N14_GLOMEINOTH",
         "E4_GRAVES_STRICT-M13_POLYARTHROPATHIES",
         "E4_GRAVES_STRICT-RHEUMA_NOS",
         "AUTOIMMUNE_NONTHYROID-K11_IBD_STRICT"
         )

ggplot(cl %>% filter(traitpair==tps[1]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[1])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[2]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[2])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[3]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[3])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic\nCorrelation","Hap Reg\nBlock 1","Hap Reg\nBlock 2","Hap Reg\nBlock 3","Phenotypic\nCorrelation"))+
  theme_classic()+theme(text=element_text(size=30))

ggplot(cl %>% filter(traitpair==tps[4]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[4])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[5]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[5])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[6]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[6])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[7]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[7])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[8]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[8])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[9]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[9])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[10]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[10])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[11]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[11])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

ggplot(cl %>% filter(traitpair==tps[12]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=tps[12])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

more <- c("M13_SJOGREN-M13_RHEUMA",
          "K11_IBD_STRICT-K11_IBD_STRICT_PSC")
ggplot(cl %>% filter(traitpair==more[1]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=more[1])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()
ggplot(cl %>% filter(traitpair==more[2]),
       aes(x=cor,y=measure,fill=measure))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(xmin=selow,xmax=sehigh),width=0.2)+ # error bars are cor +- SE
  labs(x="Correlation Value",y="Correlation Measure",title=more[2])+
  guides(fill="none")+
  scale_y_discrete(labels=c("Genetic Correlation","Hap Reg Block 1","Hap Reg Block 2","Hap Reg Block 3","Phenotypic Correlation"))+
  theme_classic()

##########################
# Chi squared
chisq.test(table(dat_bin %>% select(AB1_EBV,G6_MS)))

temp <- dat_bin %>% select(G6_MS,AUTOIMMUNE_NONTHYROID)
chisq.test(table(temp))

temp <- dat_bin %>% select(E4_GRAVES_STRICT,M13_RHEUMA)
chisq.test(table(temp))

temp <- dat_bin %>% select(E4_GRAVES_STRICT,AB1_ERYSIPELAS)
chisq.test(table(temp))

temp <- dat_bin %>% select(E4_GRAVES_STRICT,AB1_TUBERCULOSIS)
chisq.test(table(temp))


cov<-fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_PHENO_V1.txt.gz")
chisq.test(table(cov %>% select(AB1_VARICELLA,AB1_ZOSTER)))

df_full <- data.frame()

for (t1 in keep){
  for (t2 in keep) {
    if(t1==t2){next}
    c <- chisq.test(table(cov %>% select(!!sym(t1),!!sym(t2))))
    s <- c$statistic
    p <- c$p.value
    df <- data.frame(trait1=t1,trait2=t2,s=s,p=p)
    df_full=rbindlist(list(df_full,df))
  }
}

df_full <- df_full %>% mutate(traitpair=paste0(trait1,"-",trait2))

# Could make plot of trait pairs sig with chi squared but not phenocor
notsigpc <- filter(pfilt,p>0.05)$traitpair
sigchi <- filter(df_full, traitpair %in% notsigpc & p<1e-6)
ggplot(sigchi,aes(x=s,y=reorder(traitpair,s)))+
  labs(x="Chi Squared Statistic",y="Trait Pair")+
  geom_point()+theme_classic()
