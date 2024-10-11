############################### regress_enrich.R
# Courtney Smith
# Goal of script: Disease burden analysis

library(dplyr)

nsnp=1000

# set directories
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

p <- data.table::fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_PHENO_V1.txt.gz")
combstats <- list()
blocks=1:3

for (bnum in blocks){
  # For each block get the clusters and traits to keep
  r <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_regdatatoplot.txt"),header=T)
  c <- data.table::fread(paste0("hb",bnum,"_clusterdoses.txt"))

  sumstats <- list()

  # Get the traits to look at
  sigtraits <- r$V1
  pt <- p %>% select(FINNGENID,all_of(sigtraits))

  # Get the clusters to look at
  keepclust <- as.integer(colnames(r %>% select(-V1)))

  # How many people have at least one of these traits - get baseline for the block's sig traits
  nwith=sum(rowSums(pt %>% select(-FINNGENID),na.rm=TRUE)>0)
  ntotal=nrow(pt)
  frac=sum(rowSums(pt %>% select(-FINNGENID),na.rm=TRUE)>0)/nrow(pt)

  sumstats[[100]] <- data.frame(block=bnum,cluster=0,nwith=nwith,ntotal=ntotal,frac=frac,basechange="Base Rate",perofbase=100,perchangefrombase=0) # add all people as cluster 0
  basefrac=frac

  # Which people are in which clusters (at least one dose)
  for (cnum in keepclust){
    inclust <- filter_at(c,cnum+1, all_vars(.>=1))$FINNGENID
    pfilt <- pt %>% filter(FINNGENID %in% inclust) %>% select(-FINNGENID)
    nwith=sum(rowSums(pfilt,na.rm=TRUE)>0)
    ntotal=nrow(pfilt)
    frac=sum(rowSums(pfilt,na.rm=TRUE)>0)/nrow(pfilt)
    basechange=ifelse(frac>basefrac,"Above","Below")
    sumstats[[cnum]] <- data.frame(block=bnum,cluster=cnum,nwith=nwith,ntotal=ntotal,frac=frac,basechange=basechange)%>% mutate(perofbase=frac/basefrac*100,perchangefrombase=(frac-basefrac)/basefrac*100)
  }
  combstats[[bnum]] <- bind_rows(sumstats)
}

allstats <- bind_rows(combstats) %>% arrange(frac) %>% mutate(clusterblock=paste0(block,"-",cluster))

library(ggplot2)
ggplot(allstats, aes(fill=basechange))+facet_wrap(~block,scales="free_x")+
  geom_bar(stat="identity",aes(x=reorder(as.factor(as.character(clusterblock)),frac),y=frac))+
  scale_x_discrete(breaks=allstats$clusterblock,labels=allstats$cluster)+
  labs(x="Haplotype Group",y="Fraction of individuals who have at least one\nof the Block's significant traits",fill="Relative Fraction")+
  theme_classic()+theme(legend.position="top",axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5),
                        axis.title=element_text(size=16),
                        legend.text=element_text(size=20),
                        legend.title=element_text(size=18))
ggsave(paste0("/home/ivm/haplo_analysis/figures/snps",nsnp,"/diseaseburdenbar.tiff"),height=8,width=12)

# Save disease burden results; make sure privacy standards are upheld
allstats_tosave <- allstats %>% rename(hapgroup=cluster,hapgroupblockid=clusterblock) %>% select(-nwith)
allstats_tosave$ntotal <- ifelse(allstats_tosave$ntotal<50,"<50",allstats_tosave$ntotal)
data.table::fwrite(allstats_tosave, paste0(nsnp_dir,"diseaseburdenstats.txt"), row.names = F, quote = F, sep = "\t")

##### CHI SQUARED TEST
bnum=1
cq <- allstats %>% filter(block==bnum & cluster!=0) %>% mutate(nwithout=ntotal-nwith) %>% select(cluster,nwith,nwithout)
cq$cluster <- paste0("cluster",cq$cluster)
cql <- tidyr::pivot_longer(cq,cols=c(nwith,nwithout),names_to="type",values_to="value")
cqw <- as.data.frame(tidyr::pivot_wider(cql,names_from=cluster,values_from=value) %>% select(-type))
rownames(cqw) <- c("nwith","nwithout")
stats::chisq.test(cqw)

bnum=2
cq <- allstats %>% filter(block==bnum & cluster!=0) %>% mutate(nwithout=ntotal-nwith) %>% select(cluster,nwith,nwithout)
cq$cluster <- paste0("cluster",cq$cluster)
cql <- tidyr::pivot_longer(cq,cols=c(nwith,nwithout),names_to="type",values_to="value")
cqw <- as.data.frame(tidyr::pivot_wider(cql,names_from=cluster,values_from=value) %>% select(-type))
rownames(cqw) <- c("nwith","nwithout")
stats::chisq.test(cqw)

bnum=3
cq <- allstats %>% filter(block==bnum & cluster!=0) %>% mutate(nwithout=ntotal-nwith) %>% select(cluster,nwith,nwithout)
cq$cluster <- paste0("cluster",cq$cluster)
cql <- tidyr::pivot_longer(cq,cols=c(nwith,nwithout),names_to="type",values_to="value")
cqw <- as.data.frame(tidyr::pivot_wider(cql,names_from=cluster,values_from=value) %>% select(-type))
rownames(cqw) <- c("nwith","nwithout")
stats::chisq.test(cqw)

#########################
#
bnum=1
co <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_heatmapclustorder.txt"),sep=",")
co$clustnum <- as.integer(gsub("Cluster ","",co$cluster))

as <- allstats %>% filter(block==bnum&cluster!=0)
as$cluster <- factor(as$cluster,levels=rev(co$clustnum))
as$perchangefrombase <- as.numeric(ifelse(as$perchangefrombase>5,5,as$perchangefrombase))
ggplot(as, aes(x = factor(block), y = cluster, fill = perchangefrombase)) +
  #geom_tile(color="black",size=1)+
  geom_tile(width=0.2)+
  #scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
  scale_fill_gradient2(low="green",mid="white",high="orange",midpoint=0,limits=c(-5,5))+ #limits=c(-5,5)
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  #labs(x="Traits",y="Clusters",fill="Z-score")
  labs(x="",y="",fill="Percent Change from Baseline")

bnum=2
co <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_heatmapclustorder.txt"),sep=",")
co$clustnum <- as.integer(gsub("Cluster ","",co$cluster))

as <- allstats %>% filter(block==bnum&cluster!=0)
as$cluster <- factor(as$cluster,levels=co$clustnum)
ggplot(as, aes(x = factor(block), y = cluster, fill = perchangefrombase)) +
  #geom_tile(color="black",size=1)+
  geom_tile(width=0.2)+
  #scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
  scale_fill_gradient2(low="green",mid="white",high="orange",midpoint=0)+ #limits=c(-5,5)
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  #labs(x="Traits",y="Clusters",fill="Z-score")
  labs(x="",y="",fill="Percent Change from Baseline")

bnum=3
co <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_heatmapclustorder.txt"),sep=",")
co$clustnum <- as.integer(gsub("Cluster ","",co$cluster))

as <- allstats %>% filter(block==bnum&cluster!=0)
as$cluster <- factor(as$cluster,levels=co$clustnum)
ggplot(as, aes(x = factor(block), y = cluster, fill = perchangefrombase)) +
  #geom_tile(color="black",size=1)+
  geom_tile(width=0.2)+
  #scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
  scale_fill_gradient2(low="green",mid="white",high="orange",midpoint=0)+ #limits=c(-5,5)
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  #labs(x="Traits",y="Clusters",fill="Z-score")
  labs(x="",y="",fill="Percent Change from Baseline")

bnum=1
co <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_heatmapclustorder.txt"),sep=",")
co$clustnum <- as.integer(gsub("Cluster ","",co$cluster))

###
binom.test(562,770,0.675) # 562 from 0.73*770 bc 73% of people * 770 haplotypes compared to 67.5% baseline prevalence
as <- allstats %>% filter(block==bnum&cluster!=0)
as$cluster <- factor(as$cluster,levels=co$clustnum)
as$perchangefrombase <- as.numeric(ifelse(as$perchangefrombase>5,5,as$perchangefrombase))
ggplot(as, aes(x = factor(block), y = cluster, fill = perchangefrombase)) +
  #geom_tile(color="black",size=1)+
  geom_tile(width=0.2)+
  #scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
  scale_fill_gradient2(low="green",mid="white",high="orange",midpoint=0,limits=c(-5,5))+ #limits=c(-5,5)
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  #labs(x="Traits",y="Clusters",fill="Z-score")
  labs(x="",y="",fill="Percent Change from Baseline")
