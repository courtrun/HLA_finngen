# Courtney Smith - Haplotype analysis
# Defining the blocks for the haplotype analysis, listing all snps in each and plotting plots on main horizontal LD figure

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

#############
Plot boundaries

# Get min and max from the heatmap data snp positions
pos <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/matrix/filt428hits_wide.txt")  %>% tidyr::separate(ID,sep="_",into=c("CHR","POS","A1","A2")) %>% select(POS)
pos$POS <- as.numeric(pos$POS)

# reformat data
segments_data <- data.frame(
  x = c(subregion1_bounds[1], subregion2_bounds[1], subregion3_bounds[1]),
  xend = c(subregion1_bounds[2], subregion2_bounds[2], subregion3_bounds[2]),
  y = rep(0, 3),
  yend = rep(0, 3),
  color = c("#E31A1C", "#FF7F00", "#33A02C")
)

ggplot(segments_data, aes(x = x, xend = xend, y = y, yend = yend, color = color)) +
  geom_segment(size = 3) +
  scale_color_identity() +
  labs(x="Genome Position")+
  xlim(min(pos$POS), max(pos$POS)) +
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
ggsave("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/hapblockboundaries.png")
ggsave("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/hapblockboundaries.pdf")

#############


# Create a ggplot object with the original data as lines
p <- ggplot(df, aes(x = MIN_BP, y = factor(row.names(df)))) +
  geom_segment(aes(xend = MAX_BP, yend = factor(row.names(df))), size = 1, alpha = 0.7) +
  ylab("Line Number")

#brewer.pal(n = 8, name = "Dark2")
#brewer.pal(n = 8, name = "Paired")


# Add the lines for each subregion in a different color
p + geom_segment(aes(x = subregion1_bounds[1], xend = subregion1_bounds[2], y = 0, yend = 0),
                 color = "#E31A1C", size = 2) +
  geom_segment(aes(x = subregion2_bounds[1], xend = subregion2_bounds[2], y = 0, yend = 0),
               color = "#FF7F00", size = 2) +
  geom_segment(aes(x = subregion3_bounds[1], xend = subregion3_bounds[2], y = 0, yend = 0),
               color = "#33A02C", size = 2)

s <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/AB1_ACTINOMYCOSIS.gz") # get snp info, can be any trait
s <- s %>% filter(af_alt>0.01&af_alt<0.99)

snp1 <- filter(s,`#chrom`==6 & pos > subregion1_bounds[1] & pos < subregion1_bounds[2])
snp2 <- filter(s,`#chrom`==6 & pos > subregion2_bounds[1] & pos < subregion2_bounds[2])
snp3 <- filter(s,`#chrom`==6 & pos > subregion3_bounds[1] & pos < subregion3_bounds[2])
snp2b <- filter(s,`#chrom`==6 & pos > subregion2b_bounds[1] & pos < subregion2b_bounds[2])
snp3b <- filter(s,`#chrom`==6 & pos > subregion3b_bounds[1] & pos < subregion3b_bounds[2])

# List of all genes that were annotated as a nearest gene for any snps in the block
unique(snp1$nearest_genes) #
unique(snp2$nearest_genes) #
unique(snp3$nearest_genes) #
unique(snp2b$nearest_genes) #
unique(snp3b$nearest_genes) #

# List of genes that are at least partially overlapped by the block as defined by crosses any part of within the gene boudaneries
g1 <- (g %>% filter((stop>subregion1_bounds[1]&stop<subregion1_bounds[2]) | (start>subregion1_bounds[1]&start<subregion1_bounds[2])))$nearest_genes
g2 <- (g %>% filter((stop>subregion2_bounds[1]&stop<subregion2_bounds[2]) | (start>subregion2_bounds[1]&start<subregion2_bounds[2])))$nearest_genes
g3 <- (g %>% filter((stop>subregion3_bounds[1]&stop<subregion3_bounds[2]) | (start>subregion3_bounds[1]&start<subregion3_bounds[2])))$nearest_genes

# Length of block and number of snps in block
subregion1_bounds[2]-subregion1_bounds[1]+1
length(unique(snp1$rsids))

subregion2_bounds[2]-subregion2_bounds[1]+1
length(unique(snp2$rsids))

subregion3_bounds[2]-subregion3_bounds[1]+1
length(unique(snp3$rsids))

subregion2b_bounds[1]
subregion2b_bounds[2]
subregion2b_bounds[2]-subregion2b_bounds[1]+1
length(unique(snp2b$rsids))

subregion3b_bounds[1]
subregion3b_bounds[2]
subregion3b_bounds[2]-subregion3b_bounds[1]+1
length(unique(snp3b$rsids))

set.seed(123)
snp <- gsub(",.*", "", snp1$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg1_full.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp2$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg2_full.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp3$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg3_full.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp2b$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg2b_full.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp3b$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg3b_full.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

##############################

library(dplyr)
library(ggplot2)

# LD blocks for each hit
df <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_LD08.txt")
g <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/gene_limits_annotated.txt")

## Identify gene boundaries of each gene of interest
gene1="HLA-B"
subregion1_bounds <- c(filter(g,nearest_genes==gene1)$start,filter(g,nearest_genes==gene1)$stop)

gene2="HLA-DRB1"
subregion2_bounds <- c(filter(g,nearest_genes==gene2)$start,filter(g,nearest_genes==gene2)$stop)

gene3="HLA-A"
subregion3_bounds <- c(filter(g,nearest_genes==gene3)$start,filter(g,nearest_genes==gene3)$stop)

s <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/AB1_ACTINOMYCOSIS.gz") # get snp info, can be any trait
s <- s %>% filter(af_alt>0.01&af_alt<0.99) # filter to snps with MAF > 1%

snp1 <- filter(s,`#chrom`==6 & pos >= subregion1_bounds[1] & pos <= subregion1_bounds[2])
snp2 <- filter(s,`#chrom`==6 & pos >= subregion2_bounds[1] & pos <= subregion2_bounds[2])
snp3 <- filter(s,`#chrom`==6 & pos >= subregion3_bounds[1] & pos <= subregion3_bounds[2])

# Genes this covers
unique(snp1$nearest_genes) #
unique(snp2$nearest_genes) #
unique(snp3$nearest_genes) #

# Length of block and number of snps in block
subregion1_bounds[2]-subregion1_bounds[1]+1
length(unique(snp1$rsids))

subregion2_bounds[2]-subregion2_bounds[1]+1
length(unique(snp2$rsids))

subregion3_bounds[2]-subregion3_bounds[1]+1
length(unique(snp3$rsids))

set.seed(123)
snp <- gsub(",.*", "", snp1$rsids)
OUTPUT_FILE=paste0("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/",gene1,"_full.tsv")
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp2$rsids)
OUTPUT_FILE=paste0("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/",gene2,"_full.tsv")
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp3$rsids)
OUTPUT_FILE=paste0("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/",gene3,"_full.tsv")
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)
