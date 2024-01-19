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
z <- t(z)
my_function <- function(x) {return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}
z <- t(apply(z,FUN=my_function,MARGIN=1)) # apply transform to each row (1=row, 2=column), row normalize
means <- rowMeans(as.matrix(z)) # calculate the mean of each row (aka each Trait) of z
SDs <- rowSds(as.matrix(z)) # calculate the SD of each row (aka each Trait) of z
z <- sweep(z, 1, means) # subtract each row by its respective mean
z <- sweep(z, 1, SDs, FUN = '/') # and divide each row by its respective SD
colnames(z) <- full$ID
z <- na.omit(z)
z <- t(z)
z[rowMedians(as.matrix(z))<0,] <- z[rowMedians(as.matrix(z))<0,]*-1 # flip sign so rsids have median zscore across Traits thats positive
pvalues <- full %>% select(ID,matches("pval_"))
z_full <- as.data.frame(z)
z_full$ID <- rownames(z_full)
#z.long <- melt(z, id = c("Trait"))
zsnp <- z_full %>% select(-ID)

### add spacing proportional to position in genome on x axis
pos <- z_full %>% tidyr::separate(ID,sep="_",into=c("CHR","POS","A1","A2")) %>% select(POS)
total=as.numeric(max(pos$POS))-as.numeric(min(pos$POS))
pos <- pos %>% mutate(POS=as.numeric(POS),Diff = POS - lag(POS))
pos <- pos %>% mutate(spacing=Diff*3.69207665e-5*10)  # 177.8/total which is 4815718 = mm/pos, 7 inches = 177.8 mm = width pic

library(ComplexHeatmap)

# Convert the data to a matrix if it's not already in that format
zsnpt <- as.matrix(t(zsnp))

# Create a color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

# Create a Heatmap object
splits = ncol(zsnpt)
distv <- na.omit(pos$spacing)
ht <- Heatmap(zsnpt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
column_gap=unit(distv, "mm"),column_split=c(1:splits),cluster_column_slices=FALSE,
show_column_names = FALSE, show_row_names=FALSE,column_title=NULL,
column_names_gp = gpar(fill = "gray"), show_row_dend = FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_pos.png",width=6000,height=3000)
ht
dev.off()



ht <- Heatmap(zsnpt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,
column_names_gp = gpar(fontsize = 8),
show_column_names = TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 8))

# ht = draw(ht); row_order(ht)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_cht.png",width=4000,height=4000)
ht
dev.off()

ht <- Heatmap(zsnpt,col = my_palette, cluster_rows = TRUE,cluster_columns = TRUE,
column_names_gp = gpar(fontsize = 8),
show_column_names = TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 8))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_clust_cht.png",width=4000,height=4000)
ht
dev.off()

ang <- data.table::fread(/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_condit_results.txt)
head(ang %>% arrange(pos) %>% select(pos,nearest_genes,annot))
head(ang %>% filter(grepl("3263",pos)) %>% arrange(pos) %>% select(pos,nearest_genes,annot))

zsnpt <- as.data.frame(zsnpt)
zsnpt$order <- row_order(ht)
zsnpt$traits <- rownames(zsnpt)
cols <- zsnpt %>% select(order,traits) %>% arrange(order)

an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories.txt") # Trait categories

count <- left_join(t,an,by=c("trait"))
count <- count %>% arrange(Category,-n)
count$trait <- factor(count$trait, levels=unique(count$trait))

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


# Create a Heatmap object
splits = ncol(zsnpt)
distv <- na.omit(pos$spacing)
ht <- Heatmap(zsnpt,col = my_palette, cluster_rows = TRUE,cluster_columns = FALSE,column_gap=unit(distv, "mm"),column_split=c(1:splits),cluster_column_slices=FALSE,
show_column_names = TRUE, show_row_names=FALSE)
png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/figures/428fullhits_alltraits_spaced_pos.png",width=4000,height=4000)
ht
dev.off()
