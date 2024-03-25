# Courtney Smith and Satu Strausz - Main Analysis
# Make a barplot with x-axis is gene (ordered by genome position) for each gene in the HLA region
# and y-axis is the number of hits in that region, colored by trait type that the snp was a hit for
# Also make plot with vertical lines for HLA genes and horizontal for LD boundaries for scatterplot of hits

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)

r <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/external_data/41586_2015_BFnature16152_MOESM270_ESM.txt")
h <- filter(r,CHR==6&corrected.p<5e-8&POS>28510120&POS<33480577)
hits <- h %>% mutate(SNP=ifelse(ID %in% c("rs2269424","rs2517713", "rs4713462", "rs3130187"),"Top SNP","Other"))

# Load files
lr<-fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_LD08.txt") # LD regions per SNP, r2 > 0.8
##lr<-fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_LD05.txt") # LD regions per SNP, r2 > 0.5
s <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_condit_results.txt",fill=TRUE) # HLA hits
#an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories.txt") # Trait categories
an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories_v2.tsv") # Trait categories
# Munge
an$Category <- ifelse(an$Category=="Infection","Infectious",
      ifelse(an$Category=="Neoplasm","Neoplastic",ifelse(an$Category=="Neuro","Neurologic",an$Category)))

h <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits.txt") # Gene boundaries

an <- an %>% select(-LONGNAME)
h <- h %>% mutate(gene_length=stop-start+1) # Calculate gene length
sfilt <- s %>% select(-nearest_genes,-pos,-round)
h <- left_join(sfilt,h,by=c("rsids")) %>% left_join(an,by=c("trait"))
g <- as.data.frame(h %>% group_by(nearest_genes,Category) %>% count())
g <- left_join(g,distinct(h %>% select(nearest_genes,pos_mean,gene_length,start,stop)),by=c("nearest_genes"))

s2<-merge(s, lr, by.x = "pos", by.y = "BP_A", all.y = T)
s2<-s2[!duplicated(s2$rsids),]
k<-left_join(left_join(s2,an,by=c("trait")),h %>% select(rsids, nearest_genes, start, stop, pos_mean, gene_length) %>% rename(gene_start=start, gene_stop=stop,gene_pos_mean=pos_mean), by=c("rsids"))

l <- distinct(g %>% group_by(nearest_genes,pos_mean) %>% summarize(total_n=sum(n),cm=cumsum(n)))
g <- bind_cols(g,l %>% ungroup() %>% select(-nearest_genes,-pos_mean)) %>% mutate(ymin=cm-n,ymax=ymin+n)
lb <- distinct(l %>% select(nearest_genes,pos_mean,total_n))
lb$genelabel <- ifelse(lb$total_n>2,lb$nearest_genes,"")

g$class <- factor(g$Category,levels=c("Autoimmune","Cardiometabolic","Infectious","Neoplastic","Neurologic","Organ","Other","Rheumatologic")) # Order categories for plot

# Plot
p <- ggplot(g) +
  geom_rect(aes(xmin = start/1e6, xmax = stop/1e6, ymin = ymin, ymax = ymax, fill = class)) +
  xlab("Genome Position (Mb)") +
  ylab("Number Hits") +
  labs(fill="Trait Group")+
  geom_label_repel(data=lb,aes(x=pos_mean/1e6,y=total_n,label=paste0(genelabel)),size=5) + #, min.segment.length=unit(0,'lines'), show.legend = FALSE,force=10,segment.size=0.25)+
  theme_classic()+
  scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Cardiometabolic","Infectious","Neoplastic","Neurologic","Organ","Other","Rheumatologic"))+
  theme(text = element_text(size = 16))#+  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  #+ scale_y_reverse()

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_all_counts.png", width = 1200, height = 500)
p
dev.off()

k2 <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_hla_only.txt")

k<-merge(k, k2, by="pos", all.x=T)
k<-k[!duplicated(k$pos),]
k<-k[order(k$Category),]
k$P2<-1:nrow(k)

k$Category <- factor(k$Category,levels=c("Autoimmune","Cardiometabolic","Infectious","Neoplastic","Neurologic","Organ","Other","Rheumatologic"))

png("R10_LD_plot.png", width = 1200, height = 500)
ggplot(k) +
  geom_rect(aes(xmin = MIN_BP, xmax = MAX_BP, ymin = P2-0.75, ymax = P2+0.75, fill = Category)) +
  geom_rect(aes(xmin = start2, xmax = stop2, ymin = 0, ymax = 430))+
  xlab("Genome Position of SNP LD > 0.8 region") +
  ylab("") +
  scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Cardiometabolic","Infectious","Neoplastic","Neurologic","Organ","Other","Rheumatologic"))+
  scale_color_brewer(palette = "Paired")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16))+
  ## geom_point(aes(x=pos, y=P2), shape=21, size=2)+
  geom_point(aes(x=pos, y=P2, fill=Category), size=3, color="black", pch=21)
dev.off()

png("R10_LD_plot2.png", width = 1200, height = 1500)
k %>% mutate(temp= fct_reorder(rsids,MIN_BP),P2=seq(1,nrow(.)*3,3)) %>%
  ggplot() +
  geom_rect(aes(xmin = MIN_BP, xmax = MAX_BP, ymin = P2-0.75, ymax = P2+0.75, fill = Category)) +
  geom_rect(aes(xmin = start2, xmax = stop2, ymin = 0, ymax = 430*3))+
  xlab("Genome Position of SNP LD > 0.8 region") +
  ylab("") +
  scale_fill_brewer(palette = "Paired",breaks=c("Other","Organ","Neuro","Neoplasm","MSK","Cardiometabolic","Infection","Autoimmune"))+
  scale_color_brewer(palette = "Paired")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16))+
  ## geom_point(aes(x=pos, y=P2), shape=21, size=2)+
  geom_point(aes(x=pos, y=P2, fill=Category), size=3, color="black", pch=21)
dev.off()
