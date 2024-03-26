# Courtney Smith - Main analysis
# Plot genes labeled in HLA region (x-axis is genome position and y-axis is nothing,
# colored by whether the gene is classical hla, non-classical hla, or non-hla)

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggrepel)
library(RColorBrewer)
options(ggrepel.max.overlaps = Inf)

# Load files
h <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_annotated.txt") # Gene boundaries
#h$y <- runif(n=nrow(h),min=1,max=10)
h$y <- 1
h$min_stop <- ifelse(h$gene_length<14000,h$start+15000,h$stop) # make sure each is big enough to be visible
h$HLA_gene <- ifelse(h$HLA_gene=="N","Non-HLA",ifelse(h$HLA_gene=="Classical","Classical HLA","Non-Classical HLA"))
class1list <- c("HLA-A","HLA-B","HLA-C")
h$class <- ifelse(h$HLA_gene!="Classical HLA","",ifelse(h$nearest_genes %in% class1list,"HLA Class 1","HLA Class 2"))
h$color <- ifelse(h$HLA_gene=="Classical HLA",h$class,h$HLA_gene)
h$label <- ifelse(h$class=="","",h$nearest_genes)

# Plot all protein coding genes in HLA region that are the nearest gene to at least one of the hits
p <- ggplot(h) +
  geom_rect(aes(xmin = start, xmax = min_stop, ymin = -1, ymax = y, fill = color)) +
  labs(x="Genome Position",fill="Gene Classification") +
  ylim(-20,20)+
  geom_label_repel(aes(x=pos_mean,y=y,label=label,color=color),size=12,force=50, show.legend = F) + #, min.segment.length=unit(0,'lines'), show.legend = FALSE,force=10,segment.size=0.25)+
  theme_classic()+
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(),
    axis.line=element_blank(),
    axis.title.x = element_blank(),
      axis.text.x=element_blank(),
      axis.ticks=element_blank())+
  theme(text = element_text(size = 40))+
  scale_color_manual(values=c("#D95F02","orange2","lightgreen","darkgray"))+
  scale_fill_manual(values=c("#D95F02","orange2","lightgreen","darkgray"))+
  #scale_color_brewer(palette = "Dark2")+
  #scale_fill_brewer(palette = "Dark2")+
  guides(color= guide_legend(override.aes = aes(label = "")))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/gene_plot.png", width = 1800, height = 500)
p
dev.off()

#### Repeat above but plot all protein coding genes in HLA region (not just those that are nearest gene to one of the hits)
h <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/HLA_coding_genes.txt") %>%
rename(nearest_genes=name) %>% mutate(pos_mean=(start+stop)/2,gene_length=(stop-start+1))# Gene boundaries
h$HLA_gene <- ifelse(!grepl("HLA",h$nearest_genes),"N",ifelse(grepl("HLA-A|HLA-B|HLA-C|HLA-DR|HLA-DQ|HLA-DP",h$nearest_genes),"Classical","Non-Classical"))
h$y <- 1
h$min_stop <- ifelse(h$gene_length<4000,h$start+4000,h$stop) # make sure each is big enough to be visible
h$HLA_gene <- ifelse(h$HLA_gene=="N","Not HLA",ifelse(h$HLA_gene=="Classical","Classical HLA","Non-Classical HLA"))

# Merge with genes that have a hit with it as the nearest_genes
hitgenes <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_annotated.txt")$nearest_genes # Gene boundaries
h$hitgenes <- ifelse(h$nearest_genes %in% hitgenes,"Nearest Gene for a Hit","Not Hit")

# Plot all genes (regardless of whether nearest to a hit or not), coloring by if HLA or not
p <- ggplot(h) +
  geom_rect(aes(xmin = start, xmax = min_stop, ymin = -1, ymax = y, fill = HLA_gene)) +
  xlab("Genome Position") +
  ylim(-20,20)+
  geom_label_repel(aes(x=pos_mean,y=y,label=nearest_genes,color=HLA_gene),size=5,force=30, show.legend = F) + #, min.segment.length=unit(0,'lines'), show.legend = FALSE,force=10,segment.size=0.25)+
  theme_classic()+
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(text = element_text(size = 16))+scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  guides(color= guide_legend(override.aes = aes(label = "")))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/gene_plot_full.png", width = 1800, height = 500)
p
dev.off()
