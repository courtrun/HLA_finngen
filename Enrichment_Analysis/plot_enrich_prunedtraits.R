### Courtney Smith - Enrichment Analysis - Plotting
### Started on 3-5-2023
### Goal: Make a plot of enrichment and hits in HLA vs outside HLA for traits

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

options(ggrepel.max.overlaps = Inf)

### Barplot enrichment by group
e <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrichHLA_bygroup_finngenan_prunedtraits.txt")
e <- e %>% mutate(se_HLA=sd_HLA_hits/sqrt(n),se_non_HLA=sd_non_HLA_hits/sqrt(n))

e$group <- ifelse(e$group=="Rheumatologic","Rheumatic",ifelse(e$group=="Rheumatologic Comorbidities","Rheumatic Comorbidities",e$group))

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
ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup",ntraits,"_prunedtraits.png"),
width=7,height=6)

abbr <- e20
abbr$enrichabr <- ifelse(e20$enrich>100,100,e20$enrich)
ggplot(abbr,aes(y=forcats::fct_reorder(group,-order),x=enrichabr,fill=group))+
	xlim(0, max(abbr$enrichabr)) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup",ntraits,"_noinfections_prunedtraits.png"),
width=7,height=6)

# Save table with labels for each group with the plotname labels I created
OUTPUT_FILE="/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrichHLA_bygroup_annot_prunedtraits.txt"
write.table(e %>% select(group,enrich,n,HLA_hits,non_HLA_hits,se_HLA,se_non_HLA), OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

# Make plot again for the same trait groups (regardless or n traits in the trait group) that were plotted in the original analysis after filtering on those trait groups for > 15 traits
efull <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrichHLA_bygroup_annot.txt")
eog <- efull %>% filter(n>15)
ept <- filter(e,group %in% eog$group) %>% arrange(group)
ept$order <- c(1:nrow(ept))
ggplot(ept,aes(y=forcats::fct_reorder(group,-order),x=enrich,fill=group))+
	xlim(0, max(ept$enrich)+12) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup_ogtraitgroups_prunedtraits.png"),
width=7,height=6)

# Make a plot comparing enrich original vs in the plot above
ecomb <- left_join(eog,ept,by=c("group"), suffix = c(".og", ".pruned"))
correlation <- cor(ecomb$enrich.og, ecomb$enrich.pruned, use = "complete.obs")  # Pearson correlation
r2 <- correlation^2  # R^2 value
cor.test(ecomb$enrich.og, ecomb$enrich.pruned)
ggplot(ecomb, aes(x = enrich.og, y = enrich.pruned, color = group)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +  # x = y line
  labs(
    x = "Enrichment Using All Traits",
    y = "Enrichment Using Only LDSC\nGen Cor Pruned Traits",
    color = "Trait Group",
    subtitle = paste0("R² = ", round(r2, 3), ", r = ", round(correlation, 3))
  ) +
  theme_classic()+theme(text = element_text(size=20),legend.position="none")
ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup_ogtraitgroups_vs_prunedtraits.png"),
  width=7,height=6)
###########################

# By gene
e20 <- e20 %>% mutate(enrich_gene=enrich*0.48)
ggplot(e20,aes(y=forcats::fct_reorder(group,-order),x=enrich_gene,fill=group))+
	xlim(0, max(e20$enrich_gene)+12) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_pergene_prunedtraits.png",
width=10,height=6)

# By bp
e20 <- e20 %>% mutate(enrich_bp=enrich*2.4)
ggplot(e20,aes(y=forcats::fct_reorder(group,-order),x=enrich_bp,fill=group))+
	xlim(0, max(e20$enrich_bp)+12) +
	geom_bar(stat = "identity")+theme_classic()+
	labs(y="Trait Group",x="Enrichment of HLA Region",fill="Trait Group")+
	theme(text = element_text(size=20),legend.position="none",
	axis.title.y=element_blank())
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_perbp_prunedtraits.png",
width=10,height=6)

#################

### Scatterplot of enrichment vs num hits
h <- ggplot(e20,aes(x=HLA_hits,y=enrich,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="HLA hits",y="Enrichment in HLA",color="Trait Group")+
		geom_errorbar(aes(xmin=HLA_hits-se_HLA, xmax=HLA_hits+se_HLA))+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

n <- ggplot(e20,aes(x=non_HLA_hits,y=enrich,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="Genome outside HLA hits",y="Enrichment in HLA",color="Trait Group")+
		geom_errorbar(aes(xmin=non_HLA_hits-se_non_HLA, xmax=non_HLA_hits+se_non_HLA))+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

ggarrange(h,n,ncol = 2, nrow = 1)
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_hits_prunedtraits.png",height=7,width=12)
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_hits_prunedtraits.pdf",height=7,width=12)

##

### Scatterplot of HLA hits vs non HLA hits
h <- ggplot(e20,aes(x=HLA_hits,y=non_HLA_hits,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="HLA hits",y="Genome outside HLA hits",color="Trait Group")+
		geom_errorbar(aes(xmin=HLA_hits-se_HLA, xmax=HLA_hits+se_HLA,ymin=non_HLA_hits-se_non_HLA, ymax=non_HLA_hits+se_non_HLA))+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_comparehits_prunedtraits.png",height=7,width=12)
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_comparehits_prunedtraits.pdf",height=7,width=12)

### Scatterplot of HLA hits vs non HLA hits
h <- ggplot(e20,aes(x=HLA_hits/HLA_SNP,y=non_HLA_hits/genome_SNP,color=group))+geom_point()+
geom_label_repel(aes(label=group,color=group), size=4, force=10,show.legend=FALSE)+
		theme_classic()+theme(legend.position="none")+
		labs(x="Fraction HLA hits",y="Fraction genome outside HLA hits",color="Trait Group")
ggsave("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/figures/enrichHLA_bygroup20_comparefrachits_prunedtraits.png")
