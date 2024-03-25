# Courtney Smith - Main Analysis
# Started 3-17-2023
# Make a barplot with distribution of hits/traits

library(dplyr)
library(ggplot2)
library(ggridges)

l <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_condit_results.txt")

s <- l %>% group_by(rsids) %>% count() %>% arrange(-n)
summary((s)$n)

t <- l %>% group_by(trait) %>% count() %>% arrange(-n)
summary((t)$n)

an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories.txt") # Trait categories
count <- left_join(t,an,by=c("trait"))
count <- count %>% arrange(Category,-n)
count$trait <- factor(count$trait, levels=unique(count$trait))

### RIDGELINE
s <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_condit_results.txt",fill=TRUE) # HLA hits
an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories_v2.tsv") %>%
select(-Organ)# Trait categories
# Munge
an$Category <- ifelse(an$Category=="Infection","Infectious",
      ifelse(an$Category=="Neoplasm","Neoplastic",ifelse(an$Category=="Neuro","Neurologic",an$Category)))
h <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_annotated.txt")

# Annotate genes
class1 <- c("HLA-B", "HLA-A", "HLA-C")
class2 <- c(unique(h$nearest_genes[grepl("HLA-DQ",h$nearest_genes)]),unique(h$nearest_genes[grepl("HLA-DR",h$nearest_genes)]),unique(h$nearest_genes[grepl("HLA-DP",h$nearest_genes)]))

# annotate
s <- left_join(s,an,by=c("trait")) %>% select(trait,rsids,pos,nearest_genes,Category)
s$type <- ifelse(s$nearest_genes %in% class1, "class1", ifelse(s$nearest_genes %in% class2, "class2", "non-classical"))

# Create a ridgeline plot
ggplot(s, aes(x = pos, y = Category, fill= Category)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(x="Genome Position",y="")+
  theme_minimal()
ggsave("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/ridgeline_cat.png",bg="white")


# Plot
p <- ggplot(count,aes(x=Category,y=n,fill=Category)) +
  geom_violin(width=1.1, show.legend = FALSE)+
  geom_boxplot(width=0.1, show.legend = FALSE)+
  xlab("Trait Group") +
  ylab("Number of HLA hits") +
  labs(fill="Trait Group")+
  theme_classic()+
  scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Infection","Cardiometabolic","MSK","Neoplasm","Neuro", "Organ","Other"))+
  theme(text = element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/trait_hits_violin.png",width=1500)
p
dev.off()

# Plot
countfilt <- count %>% filter(n>3)
p <- ggplot(countfilt,aes(x=trait,y=n,fill=Category)) +
  geom_bar(stat="identity") +
  xlab("Traits") +
  ylab("Number of HLA hits") +
  labs(fill="Trait Group")+
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)+
  #geom_label_repel(aes(x=Category,y=n,label=paste0(n),fill="white"),size=5, show.legend = FALSE) + #, min.segment.length=unit(0,'lines'),force=10,segment.size=0.25)+
  theme_classic()+
  scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Infection","Cardiometabolic","MSK","Neoplasm","Neuro", "Organ","Other"))+
  theme(text = element_text(size = 16))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/trait_hitsfilt250.png",width=1200)
p
dev.off()

p <- ggplot(s%>%mutate(snp="HLA hits"),aes(x=snp,y=n)) +
  geom_violin()+
  geom_boxplot(width=0.1)+
  #xlab("HLA hits") +
  ylab("Number of Associated traits") +
  theme_classic()+
  theme(text = element_text(size = 25))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/hit_dist_violin.png")
p
dev.off()

###
p <- ggplot(count,aes(x=trait,y=n,fill=Category)) +
  geom_bar(stat="identity") +
  xlab("Traits") +
  ylab("Number of HLA hits") +
  labs(fill="Trait Group")+
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25, size=2)+
  #geom_label_repel(aes(x=Category,y=n,label=paste0(n),fill="white"),size=5, show.legend = FALSE) + #, min.segment.length=unit(0,'lines'),force=10,segment.size=0.25)+
  theme_classic()+
  scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Infection","Cardiometabolic","MSK","Neoplasm","Neuro", "Organ","Other"))+
  theme(text = element_text(size = 16))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) # element_text(angle = 45, vjust = 1, hjust=1))

#png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/trait_hits.png",width=1500)
#p
#dev.off()

##
ggplot(s, aes(x = pos, y = type, fill = type)) + facet_wrap(~Category)+
  geom_density_ridges(alpha = 0.7) +
  scale_fill_viridis_d() +
  theme_minimal()
ggsave("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/ridgeline_typecat.png",bg="white")
