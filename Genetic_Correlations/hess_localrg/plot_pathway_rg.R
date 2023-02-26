### Courtney Smith - Local Genetic Correlation Analysis - Plot local pathway values
### Started on 2-25-2023

library(dplyr)
library(ggplot2)

setwd("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19")
pgr <- data.table::fread("pathway_group_rhog.txt")
pgh <- data.table::fread("pathway_group_h2.txt")

pgr$p1 <- gsub("-.*", "", pgr$V1)
pgr$p2 <- gsub(".*-", "", pgr$V1)
pgr <- inner_join(pgr, pgh, by=c("p1" = "V1", "V2"), suffix=c("", ".p1"))
pgr <- inner_join(pgr, pgh, by=c("p2" = "V1", "V2"), suffix=c("", ".p2"))

# Calculate local rhog - rg = cov_g(A,B) / sqrt(h2(A)*h2(B))
pgr <- pgr %>% mutate(t=V4/sqrt(V4.p1^2 + V4.p2^2),se=V5/sqrt(V4.p1^2 + V4.p2^2)) %>% mutate(pval=2*pnorm(-abs(t/se))) # p = exp(-0.717*(t/se) - 0.416*(t/se)^2)

# Analyze results
l <- tidyr::spread(pgr %>% select(V1,V2,t), key = V2, value = t)
head(filter(l,sign(CLASS_I)!=sign(CLASS_II)) %>% arrange(-abs(CLASS_I)))
head(filter(l,sign(HLA)!=sign(GENOME)) %>% arrange(-abs(HLA)))

pgr$V2 <- ifelse(pgr$V2=="REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","REACT_I",ifelse(pgr$V2=="REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION","REACT_II",pgr$V2))

# Plot results
pdf("localrhog_plots.pdf", height=4, width=6, useDingbats=F);
for (pair in unique(pgr$V1)) {
    plot <- filter(pgr, V1 == pair, V2 %in% c("NONE", "ALL", "NONE_HLA", "GENOME", "CLASS_I", "CLASS_II", "HLA_REGION", "REACT_I", "REACT_II", "REACTOME")) %>% arrange(V4/V5) %>% mutate(V2 = factor(V2, levels=c("NONE", "ALL", "NONE_HLA", "GENOME", "CLASS_I", "CLASS_II", "HLA_REGION", "REACT_I", "REACT_II", "REACTOME"))) %>%
     ggplot(aes(x=reorder(V2,V4), y=t, ymin=t-se, ymax=t+se, fill=V2)) + geom_bar(stat="identity") + geom_errorbar() + theme_bw() + scale_fill_brewer(palette="Dark2") + coord_flip() + guides(fill=F) + ylab("Pathway")+
     xlab("Local rg") + ggtitle(paste0(pair))
    print(plot)
}
dev.off()
