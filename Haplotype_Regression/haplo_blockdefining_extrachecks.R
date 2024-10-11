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

#############
# Get min and max from the heatmap data snp positions
pos <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/matrix/filt428hits_wide.txt")  %>% tidyr::separate(ID,sep="_",into=c("CHR","POS","A1","A2")) %>% select(POS)
pos$POS <- as.numeric(pos$POS)

s <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/AB1_ACTINOMYCOSIS.gz") # get snp info, can be any trait
s <- s %>% filter((af_alt>0.01& af_alt<0.05) | (af_alt<0.99 & af_alt>0.95)) # filter to snps between MAF 1-5%

snp1 <- filter(s,`#chrom`==6 & pos > subregion1_bounds[1] & pos < subregion1_bounds[2])
snp2 <- filter(s,`#chrom`==6 & pos > subregion2_bounds[1] & pos < subregion2_bounds[2])
snp3 <- filter(s,`#chrom`==6 & pos > subregion3_bounds[1] & pos < subregion3_bounds[2])

# List of all genes that were annotated as a nearest gene for any snps in the block
unique(snp1$nearest_genes) #
unique(snp2$nearest_genes) #
unique(snp3$nearest_genes) #

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

set.seed(123)
snp <- gsub(",.*", "", snp1$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg1_full_maf01to05.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp2$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg2_full_maf01to05.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp3$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg3_full_maf01to05.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

###################################################
s <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/AB1_ACTINOMYCOSIS.gz") # get snp info, can be any trait
s <- s %>% filter((af_alt>0.05) & (af_alt<0.95)) # filter to snps >5%

snp1 <- filter(s,`#chrom`==6 & pos > subregion1_bounds[1] & pos < subregion1_bounds[2])
snp2 <- filter(s,`#chrom`==6 & pos > subregion2_bounds[1] & pos < subregion2_bounds[2])
snp3 <- filter(s,`#chrom`==6 & pos > subregion3_bounds[1] & pos < subregion3_bounds[2])

# List of all genes that were annotated as a nearest gene for any snps in the block
unique(snp1$nearest_genes) #
unique(snp2$nearest_genes) #
unique(snp3$nearest_genes) #

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

set.seed(123)
snp <- gsub(",.*", "", snp1$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg1_full_maf05.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp2$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg2_full_maf05.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

set.seed(123)
snp <- gsub(",.*", "", snp3$rsids)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/hbg3_full_maf05.tsv"
write.table(snp, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)
