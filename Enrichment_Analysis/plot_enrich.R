### Courtney Smith - Enrichment Analysis - Plotting
### Started on 3-5-2023
### Goal: Make a plot of enrichment and hits in HLA vs outside HLA for traits

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

options(ggrepel.max.overlaps = Inf)

### Barplot enrichment by group
e <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrichHLA_bygroup_finngenan.txt")
e <- e %>% mutate(se_HLA=sd_HLA_hits/sqrt(n),se_non_HLA=sd_non_HLA_hits/sqrt(n))

# Set ntraits, min num traits a group must have to be plotted
ntraits=15 # 15, 20

# Subset
e20 <- e %>% filter(n>=ntraits) %>% arrange(group)
e20$order <- c(1:nrow(e20))

# Make plot
ggplot(e20,aes(y=forcats::fct_reorder(group,-order),x=enrich,fill=group))+
	xlim(0, max(e20$enrich)+12) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup",ntraits,".png"),
width=7,height=6)

# Save table with labels for each group with the plotname labels I created
OUTPUT_FILE="/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrichHLA_bygroup_annot.txt"
write.table(e, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)
###########################

# By gene
e20 <- e20 %>% mutate(enrich_gene=enrich*0.48)
ggplot(e20,aes(y=forcats::fct_reorder(group,-order),x=enrich_gene,fill=group))+
	xlim(0, max(e20$enrich_gene)+12) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_pergene.png",
width=10,height=6)

# By bp
e20 <- e20 %>% mutate(enrich_bp=enrich*2.4)
ggplot(e20,aes(y=forcats::fct_reorder(group,-order),x=enrich_bp,fill=group))+
	xlim(0, max(e20$enrich_bp)+12) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_perbp.png",
width=10,height=6)

#################

### Scatterplot of enrichment vs num hits
h <- ggplot(e20,aes(x=HLA_hits,y=enrich,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="HLA hits",y="Enrichment in HLA",color="Trait Group")+
		geom_errorbar(aes(xmin=HLA_hits-se_HLA, xmax=HLA_hits+se_HLA))

n <- ggplot(e20,aes(x=non_HLA_hits,y=enrich,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="Genome outside HLA hits",y="Enrichment in HLA",color="Trait Group")+
		geom_errorbar(aes(xmin=non_HLA_hits-se_non_HLA, xmax=non_HLA_hits+se_non_HLA))

ggarrange(h,n,ncol = 2, nrow = 1)
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_hits.png")

##

### Scatterplot of HLA hits vs non HLA hits
h <- ggplot(e20,aes(x=HLA_hits,y=non_HLA_hits,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="HLA hits",y="Genome outside HLA hits",color="Trait Group")+
		geom_errorbar(aes(xmin=HLA_hits-se_HLA, xmax=HLA_hits+se_HLA,ymin=non_HLA_hits-se_non_HLA, ymax=non_HLA_hits+se_non_HLA))
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_comparehits.png")

### Scatterplot of HLA hits vs non HLA hits
h <- ggplot(e20,aes(x=HLA_hits/HLA_SNP,y=non_HLA_hits/genome_SNP,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="Fraction HLA hits",y="Fraction genome outside HLA hits",color="Trait Group")
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_comparefrachits.png")
