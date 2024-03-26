library(dplyr)
library(ggplot2)

nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

avr <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_regresults.txt")) %>% mutate(z=beta/se) # Load in all independent alleles at once reg
air <- data.table::fread(paste0(nsnp_dir,"alleles_indiv_regresults.txt")) %>% mutate(z=beta/se) # Load in each allele individually reg
cr1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cr2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression
cr3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults.txt")) %>% mutate(z=beta/se) # Load in cluster regression

ggplot(avr,aes(x=nimi,y=z))+
  geom_point()+
  labs(x="Traits",y="Z-score of VIF indep alleles")+
  theme_classic()

savr <- avr %>% filter(abs(z)>4.75)
sair <- air %>% filter(abs(z)>4.75)
scr1 <- cr1 %>% filter(abs(z)>4.75)
scr2 <- cr2 %>% filter(abs(z)>4.75)
scr3 <- cr3 %>% filter(abs(z)>4.75)

savrt <- unique(savr$nimi)
sairt <- unique(sair$nimi)
scrt1 <- unique(scr1$nimi)
scrt2 <- unique(scr2$nimi)
scrt3 <- unique(scr3$nimi)

length(unique(c(scrt1,scrt2,scrt3)))
clustsigtraits <- unique(c(scrt1,scrt2,scrt3))
clustsigonlytraits <- setdiff(clustsigtraits,savrt)

#s <- data.frame(trait=c(savrt,sairt,scrt1,scrt2,scrt3),reg=c(rep("Allele VIF Indep",length(savrt)),
 #          rep("Allele Indiv",length(sairt)),
  #         rep("Cluster Reg 1",length(scrt1)),
   #        rep("Cluster Reg 2",length(scrt2)),
    #       rep("Cluster Reg 3",length(scrt3))))
s <- data.frame(trait=c(savrt,scrt1,scrt2,scrt3),reg=c(rep("Allele VIF Indep",length(savrt)),
           rep("Cluster Reg 1",length(scrt1)),
           rep("Cluster Reg 2",length(scrt2)),
           rep("Cluster Reg 3",length(scrt3))))
s$trait <- factor(s$trait,levels=unique(s$trait))
## Comparing num of traits sig by each type of regression
s %>% group_by(reg) %>% count()
s <- left_join(s,s %>% group_by(trait) %>% count(),by=c("trait")) %>%
  arrange(n,reg)
length(unique(s$trait))
setdiff(unique(c(scrt1,scrt2,scrt3)),sairt)

# Plot Traits stacked bar plot across reg types
ggplot(s,aes(x=reorder(trait,n),fill=reg))+
  geom_bar()+labs(x="Traits",fill="Regression Type")+
  theme_classic()+theme(axis.text.x=element_text(angle=75,hjust=1))
# Plot Reg type vs traits
ggplot(s,aes(x=trait,y=reg,fill=reg))+
  geom_tile()+labs(x="Traits",fill="Regression Type")+
  theme_classic()+theme(axis.text.x=element_text(angle=75,hjust=1))
# Plot just the traits sig in at least one cluster; reg type vs traits
ggplot(s %>% filter(trait %in% clustsigtraits),aes(x=trait,y=reg,fill=reg))+
  geom_tile()+labs(x="Traits",fill="Regression Type")+
  theme_classic()+theme(axis.text.x=element_text(angle=75,hjust=1))
# Plot ONLY the traits sig in at least one cluster and not allele VIF; reg type vs traits
ggplot(s %>% filter(trait %in% clustsigonlytraits),aes(x=trait,y=reg,fill=reg))+
  geom_tile()+labs(x="Traits",fill="Regression Type")+
  theme_classic()+theme(axis.text.x=element_text(angle=75,hjust=1))

# Comparing num traits sig for a given "locus"
scr1 %>% group_by(cluster) %>% count()
sair %>% group_by(allele) %>% count()

pair <- air
w <- tidyr::spread(pair %>% select(allele,nimi,z),allele,z)
rownames(w) <- w$nimi
w <- w %>% select(-nimi)
hc <- hclust(dist(w,method="euclidean"),method="ward.D")
w <- w[hc$order,]
w$nimi <- rownames(w)
pair <- reshape2::melt(w) %>% rename(allele=variable,z=value)
pair$z <- ifelse(pair$z>5,5,ifelse(pair$z<(-5),-5,pair$z))
pair$nimi <- factor(pair$nimi,levels=unique(pair$nimi))
ggplot(pair,aes(x=allele,y=nimi,fill=z))+geom_tile()+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits = c(-5, 5))+
  theme_classic()+
  labs(x="Alleles (From Regression with one allele at a time)",y="Traits")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1),axis.text.y = element_blank())

pavr <- avr
w <- tidyr::spread(pavr %>% select(allele,nimi,z),allele,z)
rownames(w) <- w$nimi
w <- w %>% select(-nimi)
hc <- hclust(dist(w,method="euclidean"),method="ward.D")
w <- w[hc$order,]
w$nimi <- rownames(w)
pavr <- reshape2::melt(w) %>% rename(allele=variable,z=value)
pavr$z <- ifelse(pavr$z>5,5,ifelse(pavr$z<(-5),-5,pavr$z))
pavr$nimi <- factor(pavr$nimi,levels=unique(pavr$nimi))
ggplot(pavr,aes(x=allele,y=nimi,fill=z))+geom_tile()+
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits = c(-5, 5))+
  theme_classic()+
  labs(x="Alleles (From Regression with all `Independent` alleles)",y="Traits")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1),axis.text.y = element_blank())

