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
