library(dplyr)
library(ggplot2)
library(ggrepel)

nsnp="1000"
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

#### BLOCK 3 below
bnum=3

g <- data.table::fread("/home/ivm/from_satu/general_files/HLA_coding_genes.tsv")
s <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_snps.txt"))
splt <- tidyr::separate(s,V1,sep="_",into=c("CHR","POS","A1","A2"))

gene_list <- list()

for (i in 1:nrow(splt)) {
  pos <- as.numeric(splt$POS[i])
  gene <- filter(g,start <= pos & stop >= pos)$name
  gene_list[[i]] <- gene
}

splt$name <- gene_list
#splt$name <- unlist(lapply(splt$name, function(x) ifelse(length(x) == 0, NA, x)))

splt <- as.data.frame(splt %>% arrange(POS))

# Fix the issue
filter(splt,grepl(",",name))
splt$name[splt$POS==filter(splt,grepl(",",name))$POS] <- "GPSM3"
filter(splt,grepl(",",name)) # should be empty

# Add in other genes in block that didnt have snp within them
bound_lower=32094910
bound_upper=32847125
allg <- g %>% filter(start<bound_upper&stop>bound_lower) %>% arrange(start) %>%
  filter(!(name %in% splt$name))# filter to missing genes
for (gene in unique(allg$name)){# find snp immediately before the gene and add gene
  gstart <- filter(allg,name==gene)$start 
  gpos <- tail(filter(splt,POS<gstart),1)$POS
  splt[splt$POS == gpos,]$name <- gene
}

splt$order <- 1:nrow(splt)
#splt2 <- tidyr::separate_rows(splt,name,sep=", ")

filt <- as.data.frame(filter(splt,name!="character(0)") %>% # remove snps outside gene boundaries
  group_by(name) %>%
  filter(order == max(order) | order == min(order)))
gene_one <- (as.data.frame(filt %>% group_by(name) %>% count()) %>% filter(n==1))$name # filter to all genes that only have one row
filt <- bind_rows(filt,filter(filt,name %in% gene_one) %>% mutate(order=order+0.7)) %>% arrange(POS) # duplicate all rows that only have one row per gene
as.data.frame(filt %>% group_by(name) %>% count()) %>% filter(n!=2) # make sure this is empty since should be 2 for each
filt$endpoints <- rep(c("start","stop"),times=nrow(filt)/2)
filt <- filt %>% select(name,order,endpoints)
gene_pos <- tidyr::spread(as.data.frame(filt),key=endpoints,value=order)
gene_pos <- gene_pos %>% mutate(mid=(start+stop)/2)

options(ggrepel.max.overlaps = Inf)
set.seed(20)
# Create plot with all segments on single horizontal line
ggplot(gene_pos, aes(x=start, xend=stop, y=1, yend=1)) +
  geom_segment(size=5, color="darkgray") +
  labs(x = "Position", y = "") +
  geom_text_repel(aes(label = name, x=mid),y = 1, box.padding = 1, min.segment.length = unit(0, 'lines'), size = 10) + #  geom_text(aes(label = name), y = 1, vjust = -1, size = 3) +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#### BLOCK 2 below

bnum=2

g <- data.table::fread("/home/ivm/from_satu/general_files/HLA_coding_genes.tsv")
s <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_snps.txt"))
splt <- tidyr::separate(s,V1,sep="_",into=c("CHR","POS","A1","A2"))

gene_list <- list()

for (i in 1:nrow(splt)) {
  pos <- as.numeric(splt$POS[i])
  gene <- filter(g,start <= pos & stop >= pos)$name
  gene_list[[i]] <- gene
}

splt$name <- gene_list
#splt$name <- unlist(lapply(splt$name, function(x) ifelse(length(x) == 0, NA, x)))

splt <- as.data.frame(splt %>% arrange(POS))

# Fix the issue
filter(splt,grepl(",",name)) # should be empty

# Add in other genes in block that didnt have snp within them
bound_lower=31168798
bound_upper=31611071
allg <- g %>% filter(start<bound_upper&stop>bound_lower) %>% arrange(start) %>%
  filter(!(name %in% splt$name))# filter to missing genes
for (gene in unique(allg$name)){# find snp immediately before the gene and add gene
  gstart <- filter(allg,name==gene)$start 
  gpos <- tail(filter(splt,POS<gstart),1)$POS
  splt[splt$POS == gpos,]$name <- gene
}

splt$order <- 1:nrow(splt)
#splt2 <- tidyr::separate_rows(splt,name,sep=", ")

filt <- as.data.frame(filter(splt,name!="character(0)") %>% # remove snps outside gene boundaries
                        group_by(name) %>%
                        filter(order == max(order) | order == min(order)))
gene_one <- (as.data.frame(filt %>% group_by(name) %>% count()) %>% filter(n==1))$name # filter to all genes that only have one row
filt <- bind_rows(filt,filter(filt,name %in% gene_one) %>% mutate(order=order+0.7)) %>% arrange(POS) # duplicate all rows that only have one row per gene
as.data.frame(filt %>% group_by(name) %>% count()) %>% filter(n!=2) # make sure this is empty since should be 2 for each
filt$endpoints <- rep(c("start","stop"),times=nrow(filt)/2)
filt <- filt %>% select(name,order,endpoints)
gene_pos <- tidyr::spread(as.data.frame(filt),key=endpoints,value=order)
gene_pos <- gene_pos %>% mutate(mid=(start+stop)/2)

options(ggrepel.max.overlaps = Inf)
set.seed(20)

# Create plot with all segments on single horizontal line
ggplot(gene_pos, aes(x=start, xend=stop, y=1, yend=1)) +
  geom_segment(size=5, color="darkgray") +
  labs(x = "Position", y = "") +
  geom_text_repel(aes(label = name, x=mid),y = 1, box.padding = 1, min.segment.length = unit(0, 'lines'), size = 10) + #  geom_text(aes(label = name), y = 1, vjust = -1, size = 3) +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

##### BLOCK 1 below

bnum=1

g <- data.table::fread("/home/ivm/from_satu/general_files/HLA_coding_genes.tsv")
s <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_snps.txt"))
splt <- tidyr::separate(s,V1,sep="_",into=c("CHR","POS","A1","A2"))

gene_list <- list()

for (i in 1:nrow(splt)) {
  pos <- as.numeric(splt$POS[i])
  gene <- filter(g,start <= pos & stop >= pos)$name
  gene_list[[i]] <- gene
}

splt$name <- gene_list
#splt$name <- unlist(lapply(splt$name, function(x) ifelse(length(x) == 0, NA, x)))

splt <- as.data.frame(splt %>% arrange(POS))

# Fix the issue
filter(splt,grepl(",",name)) # should be empty

# Add in other genes in block that didnt have snp within them
bound_lower=29622820
bound_upper=30045616
allg <- g %>% filter(start<bound_upper&stop>bound_lower) %>% arrange(start) %>%
  filter(!(name %in% splt$name))# filter to missing genes
for (gene in unique(allg$name)){# find snp immediately before the gene and add gene
  gstart <- filter(allg,name==gene)$start 
  gpos <- tail(filter(splt,POS<gstart),1)$POS
  splt[splt$POS == gpos,]$name <- gene
}

splt$order <- 1:nrow(splt)
#splt2 <- tidyr::separate_rows(splt,name,sep=", ")

filt <- as.data.frame(filter(splt,name!="character(0)") %>% # remove snps outside gene boundaries
                        group_by(name) %>%
                        filter(order == max(order) | order == min(order)))
gene_one <- (as.data.frame(filt %>% group_by(name) %>% count()) %>% filter(n==1))$name # filter to all genes that only have one row
filt <- bind_rows(filt,filter(filt,name %in% gene_one) %>% mutate(order=order+0.7)) %>% arrange(POS) # duplicate all rows that only have one row per gene
as.data.frame(filt %>% group_by(name) %>% count()) %>% filter(n!=2) # make sure this is empty since should be 2 for each
filt$endpoints <- rep(c("start","stop"),times=nrow(filt)/2)
filt <- filt %>% select(name,order,endpoints)
gene_pos <- tidyr::spread(as.data.frame(filt),key=endpoints,value=order)
gene_pos <- gene_pos %>% mutate(mid=(start+stop)/2)

options(ggrepel.max.overlaps = Inf)
set.seed(20)
# Create plot with all segments on single horizontal line
ggplot(gene_pos, aes(x=start, xend=stop, y=1, yend=1)) +
  geom_segment(size=5, color="darkgray") +
  labs(x = "Position", y = "") +
  geom_text_repel(aes(label = name, x=mid),y = 1, box.padding = 1, min.segment.length = unit(0, 'lines'), size = 10) + #  geom_text(aes(label = name), y = 1, vjust = -1, size = 3) +
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

