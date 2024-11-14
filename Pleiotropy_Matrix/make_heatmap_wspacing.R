# Courtney Smith - HLA Project - HLA Heatmap Plot
# Script started 1-18-2023
# Goal of script: Plot heatmap of zscores for traits vs variants ascertained in them; with and without x axis

library("ggplot2")
library("ggdendro")
library("reshape2")
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggbiplot)
library(forcats)
library("pheatmap")

full <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/matrix/filt428hits_wide.txt") #

betas_full <- full %>% select(matches("^beta_"))
colnames(betas_full) <- gsub("beta_","",colnames(betas_full))
se_full <- full %>% select(matches("sebeta_"))
z <- betas_full/se_full

# Munge: (sign)*z^2/(Z max for each trait across all snps)^2
sq <- z^2
maxes <- summarise_each(abs(z), funs(max(., na.rm=TRUE)))
z <- mapply("/", sq, (maxes)^2)

# Set sign
z <- z * sign(betas_full)
z[rowMedians(as.matrix(z))<0,] <- z[rowMedians(as.matrix(z))<0,]*-1 # flip sign so rsids have median zscore across Traits thats positive

#
rownames(z) <- full$ID
z_full <- as.data.frame(z)
z_full$ID <- full$ID

### add spacing proportional to position in genome on x axis
pos <- z_full %>% tidyr::separate(ID,sep="_",into=c("CHR","POS","A1","A2")) %>% select(POS)
total=as.numeric(max(pos$POS))-as.numeric(min(pos$POS))
pos <- pos %>% mutate(POS=as.numeric(POS),Diff = POS - lag(POS))
pos <- pos %>% mutate(spacing=Diff*3.69207665e-5*10)  # 177.8/total which is 4815718 = mm/pos, 7 inches = 177.8 mm = width pic

library(ComplexHeatmap)

# Convert the data to a matrix if it's not already in that format
zt <- as.matrix(t(z))

# Create a color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

# Create a Heatmap object
splits = ncol(zt)
distv <- na.omit(pos$spacing)

# Heatmap ordered by position for snps, clustered traits, no split
ht <- Heatmap(zt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
show_column_names = TRUE, show_row_names=TRUE,column_title=NULL,
column_names_gp = gpar(fontsize=24),
row_names_gp = gpar(fontsize=24),
show_row_dend = FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_z2.png",width=6000,height=3000)
ht
dev.off()

# Same but spread by position
ht <- Heatmap(zt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
column_gap=unit(distv, "mm"),column_split=c(1:splits),cluster_column_slices=FALSE,
show_column_names = FALSE, show_row_names=FALSE,column_title=NULL,
column_names_gp = gpar(fill = "gray"), show_row_dend = FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_pos_z2.png",width=6000,height=3000)
ht
dev.off()

##########
############### Add labels for plotting - to help annotation of 3b
z_full$pos <- pos$POS
z_full <- z_full %>% arrange(pos)

# For key genes, look up gene boundaries and which SNPs mark the start and stop so can label in affinity design
g <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_annotated.txt") %>% arrange(start)
key_genes <- c("HLA-F","HLA-G","HLA-A","HLA-C","HLA-B","MICA","NOTCH4") # + "HLA-DR","HLA-DQ" genes

# Initialize an empty dataframe to store results
snp_labels <- data.frame(lower = numeric(), upper = numeric(), stringsAsFactors = FALSE)
g <- filter(g,nearest_genes %in% key_genes | grepl("HLA-DR",nearest_genes)|grepl("HLA-DQ",nearest_genes))
# For each gene in genes of interest
for (i in 1:nrow(g)) {
# Identify the start and stop positions
  start_val <- g$start[i]
  stop_val <- g$stop[i]

  # Find the lowest SNP position in z_full that is greater than start_val
  low_snp <- (z_full %>% filter(pos > start_val) %>% slice_min(pos, n = 1))$ID

  # Find the highest SNP position in z_full that is less than stop_val
  high_snp <- (z_full %>% filter(pos < stop_val) %>% slice_max(pos, n = 1))$ID

  # Append results to the results dataframe
  snp_labels <- rbind(snp_labels, data.frame(gene = g$nearest_genes[i], snp = ifelse(length(low_snp) > 0, low_snp, NA)))
  snp_labels <- rbind(snp_labels, data.frame(gene = g$nearest_genes[i], snp = ifelse(length(high_snp) > 0, high_snp, NA)))
  }

# Make labels for x axis
snp_labels %>% group_by(snp) %>% count() %>% filter(freq>1) # remove do to bug but add back manually
snp_labels %>% arrange(snp) # see which ones have repeat snps but dif genes
snp_labels <- filter(snp_labels,!(gene=="HLA-DQA2" & snp=="chr6_32763163_C_A"))
x_labels <- left_join(z_full %>% select(ID,pos) %>% distinct(),distinct(snp_labels),by=c("ID"="snp"))$gene
x_labels <- sapply(x_labels, function(x) ifelse(is.na(x), "-", x))
colnames(zt) <- x_labels

# Make y axis labels - label diseases sig in hap analysis (fig 5)
p <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/phenotype_info/new_plot_names.txt",header=T) # load in the list of traits and filter to those that are plotted
p2 <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/phenotype_info/regdatafull_long.txt",header=T) # load in the list of traits and filter to those that are plotted
custom_traits <- unique(filter(p2,plotted=="yes")$traits)
p$cols <- ifelse(p$trait %in% custom_traits,p$Plot,"-")
rownames(zt) <- p$cols

# Spread by position
ht <- Heatmap(zt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
column_gap=unit(distv, "mm"),column_split=c(1:splits),cluster_column_slices=FALSE,
show_column_names = TRUE, show_row_names=TRUE,column_title=NULL,
column_names_gp = gpar(fill = "gray"), show_row_dend = FALSE)

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_z2_spread_wlabels.png",width=6000,height=3000)
ht
dev.off()
pdf("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_z2_spread_wlabels.pdf",width=17,height=20)
ht
dev.off()

################ Repeat heatmap, unsplit, traits from fig 5 labels, SNPs annotated
colnames(zt) <- z_full$ID

ht <- Heatmap(zt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
show_column_names = TRUE, show_row_names=TRUE,column_title=NULL,
column_names_gp = gpar(fontsize=20),
row_names_gp = gpar(fontsize=20),
show_row_dend = FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_z2_wlabels.png",width=6000,height=3000)
ht
dev.off()
pdf("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_z2_wlabels.pdf",width=40,height=20)
ht
dev.off()
#################### Repeat heatmap, ordered by position, split by block

library(dplyr)
library(ggplot2)

# LD blocks for each hit
df <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_LD08.txt")
g <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_annotated.txt") %>% arrange(start)

### Divide into 3
# 100kb below HLA-F through HLA-G and to 100kb past HLA-A
subregion1_bounds <- c(filter(g,nearest_genes=="HLA-F")$start-100000,filter(g,nearest_genes=="HLA-A")$stop+100000)

# 100kb below HLA-C through HLA-B and MICA to 100kb past MICB
subregion2_bounds <- c(filter(g,nearest_genes=="HLA-C")$start-100000,filter(g,nearest_genes=="MICB")$stop+100000)

# 100kb below NOTCH4 through HLA-DRX and HLA-DQX to 100kb past HLA-DQA2
subregion3_bounds <- c(filter(g,nearest_genes=="NOTCH4")$start-100000,filter(g,nearest_genes=="HLA-DQA2")$stop+100000)

# 100kb below HLA-C through 100kb past HLA-B
subregion2b_bounds <- c(filter(g,nearest_genes=="HLA-C")$start-100000,filter(g,nearest_genes=="HLA-B")$stop+100000)

# 100kb below NOTCH4 through HLA-DRX and HLA-DQX to 100kb past HLA-DQA2
subregion3b_bounds <- c(filter(g,nearest_genes=="HLA-DRA")$start-100000,filter(g,nearest_genes=="HLA-DQB2")$stop+100000)

##

z_full$pos <- pos$POS
z_full <- z_full %>% arrange(pos)

# Add in the breaks before and after each block
bks <- c(min(which(z_full$pos<subregion1_bounds[2]&pos>subregion1_bounds[1])),
  max(which(z_full$pos<subregion1_bounds[2]&pos>subregion1_bounds[1])),
  min(which(z_full$pos<subregion2_bounds[2]&pos>subregion2_bounds[1])),
  max(which(z_full$pos<subregion2_bounds[2]&pos>subregion2_bounds[1])),
  min(which(z_full$pos<subregion3_bounds[2]&pos>subregion3_bounds[1])),
  max(which(z_full$pos<subregion3_bounds[2]&pos>subregion3_bounds[1])))
bks <- c(bks[1],diff(c(bks,nrow(z_full))))
bks <- rep(c("a", "b", "c", "d", "e", "f", "g"), times = bks)

#
ht <- Heatmap(zt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
column_split=bks,cluster_column_slices=FALSE,
show_column_names = FALSE, show_row_names=FALSE,column_title=NULL,
show_row_dend = FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_z2_brokenbyblock.png",width=6000,height=3000)
ht
dev.off()

############### Add labels for plotting

# For key genes, look up gene boundaries and which SNPs mark the start and stop so can label in affinity design
key_genes <- c("HLA-F","HLA-G","HLA-A","HLA-C","HLA-B","MICA","NOTCH4") # + "HLA-DR","HLA-DQ" genes

# Initialize an empty dataframe to store results
snp_labels <- data.frame(lower = numeric(), upper = numeric(), stringsAsFactors = FALSE)
g <- filter(g,nearest_genes %in% key_genes | grepl("HLA-DR",nearest_genes)|grepl("HLA-DQ",nearest_genes))
# For each gene in genes of interest
for (i in 1:nrow(g)) {
# Identify the start and stop positions
  start_val <- g$start[i]
  stop_val <- g$stop[i]

  # Find the lowest SNP position in z_full that is greater than start_val
  low_snp <- (z_full %>% filter(pos > start_val) %>% slice_min(pos, n = 1))$ID

  # Find the highest SNP position in z_full that is less than stop_val
  high_snp <- (z_full %>% filter(pos < stop_val) %>% slice_max(pos, n = 1))$ID

  # Append results to the results dataframe
  snp_labels <- rbind(snp_labels, data.frame(gene = g$nearest_genes[i], snp = ifelse(length(low_snp) > 0, low_snp, NA)))
  snp_labels <- rbind(snp_labels, data.frame(gene = g$nearest_genes[i], snp = ifelse(length(high_snp) > 0, high_snp, NA)))
  }

# Make labels for x axis
snp_labels %>% group_by(snp) %>% count() %>% filter(freq>1) # remove do to bug but add back manually
snp_labels %>% arrange(snp) # see which ones have repeat snps but dif genes
snp_labels <- filter(snp_labels,!(gene=="HLA-DQA2" & snp=="chr6_32763163_C_A"))
x_labels <- left_join(z_full %>% select(ID,pos) %>% distinct(),distinct(snp_labels),by=c("ID"="snp"))$gene
x_labels <- sapply(x_labels, function(x) ifelse(is.na(x), "-", x))
colnames(zt) <- x_labels

# Make y axis labels - label diseases sig in hap analysis (fig 5)
p <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/phenotype_info/new_plot_names.txt",header=T) # load in the list of traits and filter to those that are plotted
p2 <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/phenotype_info/regdatafull_long.txt",header=T) # load in the list of traits and filter to those that are plotted
custom_traits <- unique(filter(p2,plotted=="yes")$traits)
p$cols <- ifelse(p$trait %in% custom_traits,p$Plot,"-")
rownames(zt) <- p$cols

ht <- Heatmap(zt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
column_split=bks,cluster_column_slices=FALSE,
show_column_names = TRUE, show_row_names=TRUE,column_title=NULL,
 row_names_gp = gpar(fontsize = 24),
column_names_gp = gpar(fontsize=24),
show_row_dend = FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_z2_brokenbyblock_wlabels.png",width=6000,height=3000)
ht
dev.off()
pdf("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_z2_brokenbyblock_wlabels.pdf",width=17,height=20)
ht
dev.off()
