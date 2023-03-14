### Courtney Smith - Local Genetic Correlation Analysis - Plot local pathway values
### Started on 2-25-2023

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(corrplot)

setwd("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19")
pgr <- data.table::fread("pathway_group_rhog.txt")
pgh <- data.table::fread("pathway_group_h2.txt")

pgr$p1 <- gsub("-.*", "", pgr$V1)
pgr$p2 <- gsub(".*-", "", pgr$V1)
pgr <- inner_join(pgr, pgh, by=c("p1" = "V1", "V2"), suffix=c("", ".p1"))
pgr <- inner_join(pgr, pgh, by=c("p2" = "V1", "V2"), suffix=c("", ".p2"))

# Calculate local rhog - rg = cov_g(A,B) / sqrt(h2(A)*h2(B)); se = se(cov_g(A,B)) / sqrt(h2(A)*h2(B))
pgr <- pgr %>% mutate(t=V4/sqrt(abs(V4.p1*V4.p2)), se=V5/sqrt(abs(V4.p1*V4.p2))) %>% mutate(pval=2*pnorm(-abs(t/se))) # p = exp(-0.717*(t/se) - 0.416*(t/se)^2)

pgr$V2 <- ifelse(pgr$V2=="REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","REACT_I",ifelse(pgr$V2=="REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION","REACT_II",pgr$V2))
pgr$t <- ifelse(pgr$t>1,1,ifelse(pgr$t<(-1),-1,pgr$t))

# Plot results
pdf("figures/localrhog_plots.pdf", height=4, width=6, useDingbats=F);
for (pair in unique(pgr$V1)) {
    plot <- filter(pgr, V1 == pair, V2 %in% c("NONE", "ALL", "NONE_HLA", "GENOME", "CLASS_I", "CLASS_II", "HLA_REGION", "REACT_I", "REACT_II", "REACTOME")) %>% arrange(V4/V5) %>% mutate(V2 = factor(V2, levels=c("NONE", "ALL", "NONE_HLA", "GENOME", "CLASS_I", "CLASS_II", "HLA_REGION", "REACT_I", "REACT_II", "REACTOME"))) %>%
     ggplot(aes(x=reorder(V2,V2), y=t, ymin=t-se, ymax=t+se, fill=V2)) + geom_bar(stat="identity") + geom_errorbar() + theme_bw() + coord_flip() + guides(fill=F) + ylab("Pathway")+ #  scale_fill_brewer(palette="Dark2") +
     xlab("Local rg") + ggtitle(paste0(pair))
    print(plot)
}
dev.off()

# Plot results
pdf("figures/localrhog_plots_noSE.pdf", height=4, width=6, useDingbats=F);
for (pair in unique(pgr$V1)) {
    plot <- filter(pgr, V1 == pair, V2 %in% c("NONE", "ALL", "NONE_HLA", "GENOME", "CLASS_I", "CLASS_II", "HLA_REGION", "REACT_I", "REACT_II", "REACTOME")) %>% arrange(V4/V5) %>% mutate(V2 = factor(V2, levels=c("NONE", "ALL", "NONE_HLA", "GENOME", "CLASS_I", "CLASS_II", "HLA_REGION", "REACT_I", "REACT_II", "REACTOME"))) %>%
     ggplot(aes(x=reorder(V2,V2), y=t, fill=V2)) + geom_bar(stat="identity") + theme_bw() +  coord_flip() + guides(fill=F) + ylab("Pathway")+
     xlab("Local rg") + ggtitle(paste0(pair))
    print(plot)
}
dev.off()

### Choosing interesting ones to visualize
l <- tidyr::spread(pgr %>% select(V1,V2,t), key = V2, value = t)
filter(l,V1 %in% ((filter(pgr,pval<0.05)) %>% group_by(V1) %>% count() %>% filter(n>=4))$V1&sign(HLA_REGION)!=sign(NONE_HLA)) # head(filter(l,sign(HLA)!=sign(GENOME)) %>% arrange(-abs(HLA)))
head(filter(l,sign(CLASS_I)!=sign(CLASS_II)) %>% arrange(-abs(CLASS_I))) # filter(l,V1 %in% ((filter(pgr,pval<0.05)) %>% group_by(V1) %>% count() %>% filter(n>=1))$V1&sign(CLASS_I)!=sign(CLASS_II))

###
for (pathway in c("HLA","CLASS_I","CLASS_II","NONE_HLA","GENOME")){
s <- filter(pgr,V2==pathway) #
temp <- distinct(bind_rows(s %>% select(p1,p2,t),s %>% select(p2,p1,t)%>%rename(p1=p2,p2=p1)))
temp$p1 <- as.character(temp$p1)
temp$p2 <- as.character(temp$p2)
m <- tidyr::spread(temp,key=p1,value=t,fill=1) %>% select(-p2)
m <- as.matrix(m)
rownames(m) <- colnames(m)
#Heatmap(m,show_row_names=T,show_column_names=T)
png(paste("localrg_corrplot_",pathway,".png"))
corrplot(m,method='color',type="lower",tl.cex=0.75,tl.col="black",col=colorRampPalette(c("blue","white","red"))(200)) # ,order="hclust",,title=pathway
dev.off()
}
