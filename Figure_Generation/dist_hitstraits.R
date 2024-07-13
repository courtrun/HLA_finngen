# Courtney Smith - Main Analysis
# Started 3-17-2023
# Make a ridgeline plot with distribution of hits by trait group

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
  theme_minimal()+theme(text = element_text(size = 20))
ggsave("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/ridgeline_cat.png",bg="white",height=7,width=12)
