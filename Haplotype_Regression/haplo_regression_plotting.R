library(dplyr)
library(pheatmap)
library(RColorBrewer)

dir="/home/ivm/hla_files/results_100/figures/"
setwd(dir)

## BLOCK 1
h1 <- data.table::fread("block1_9clusters.txt") %>% rename(trait=nimi)
filter(h1,p<0.005)

# how many traits have at least one cluster w/ p < 0.005
length(unique(filter(h1,p<0.005)$trait)) # 70 (of the 269 tested)

# filter to a small number of traits of interest (have >= 3 clusters w/ p<0.005)
head(as.data.frame(filter(h1,p<0.005) %>% group_by(trait) %>% count() %>% arrange(-n)),15)
trait3 <- (as.data.frame(filter(h1,p<0.005) %>% group_by(trait) %>% count() %>% arrange(-n)) %>% filter(n>=3))$trait
h <- filter(h1,trait %in% trait3)
h <- h %>% mutate(z=beta/se)

beta <- tidyr::spread(h %>% select(trait,cluster,beta),cluster,beta)

z <- tidyr::spread(h %>% select(trait,cluster,z),cluster,z)
traits <- z$trait
rownames(z) <- traits
z <- z %>% select(-trait)
png(paste0(dir,"b1_traits3_nocutoff.png"))
pheatmap(z)
dev.off()
breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b1_traits3_cutoff5.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()

## add in dropped cluster w/ rescale
z <- z %>% mutate(dropped=0)
# Function to rescale row values to have mean 0
rescale_mean_zero <- function(row) {
  mean_zero <- row - mean(row)
  return(mean_zero)
}
# Apply the rescaling function to each row of the dataframe
z_full <- t(apply(z, 1, rescale_mean_zero))
rownames(z_full) <- traits
png(paste0(dir,"b1_traits3_nocutoff_full.png"))
pheatmap(z_full)
dev.off()

breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b1_traits3_cutoff5_full.png"))
pheatmap(z_full,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()

# same thing on all of traits
h1 <- h1 %>% mutate(z=beta/se)
z <- tidyr::spread(h1 %>% select(trait,cluster,z),cluster,z)
traits <- z$trait
rownames(z) <- traits
z <- z %>% select(-trait)
png(paste0(dir,"b1_alltraits_nocutoff.png"))
pheatmap(z,show_rownames=F)
dev.off()
breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b1_alltraits_cutoff5.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,show_rownames=F)
dev.off()
## add in dropped cluster w/ rescale
z <- z %>% mutate(dropped=0)
# Function to rescale row values to have mean 0
rescale_mean_zero <- function(row) {
  mean_zero <- row - mean(row)
  return(mean_zero)
}
# Apply the rescaling function to each row of the dataframe
z_full <- t(apply(z, 1, rescale_mean_zero))
rownames(z_full) <- traits
png(paste0(dir,"b1_alltraits_nocutoff_full.png"))
pheatmap(z_full,show_rownames=F)
dev.off()

breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b1_alltraits_cutoff5_full.png"))
pheatmap(z_full,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,show_rownames=F)
dev.off()

####
## BLOCK 2
h1 <- data.table::fread("block2_9clusters.txt") %>% rename(trait=nimi)
filter(h1,p<0.005)

# how many traits have at least one cluster w/ p < 0.005
length(unique(filter(h1,p<0.005)$trait)) # 120 (of the 269 tested)

# filter to a small number of traits of interest (have >= 7 clusters w/ p<0.005)
head(as.data.frame(filter(h1,p<0.005) %>% group_by(trait) %>% count() %>% arrange(-n)),15)
trait3 <- (as.data.frame(filter(h1,p<0.005) %>% group_by(trait) %>% count() %>% arrange(-n)) %>% filter(n>=7))$trait
h <- filter(h1,trait %in% trait3)
h <- h %>% mutate(z=beta/se)

beta <- tidyr::spread(h %>% select(trait,cluster,beta),cluster,beta)

z <- tidyr::spread(h %>% select(trait,cluster,z),cluster,z)
traits <- z$trait
rownames(z) <- traits
z <- z %>% select(-trait)
png(paste0(dir,"b2_traits7_nocutoff.png"))
pheatmap(z)
dev.off()
breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b2_traits7_cutoff5.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()
breaksList = seq(-10, 10, by = 1)
png(paste0(dir,"b2_traits7_cutoff10.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()

## add in dropped cluster w/ rescale
z <- z %>% mutate(dropped=0)
# Function to rescale row values to have mean 0
rescale_mean_zero <- function(row) {
  mean_zero <- row - mean(row)
  return(mean_zero)
}
# Apply the rescaling function to each row of the dataframe
z_full <- t(apply(z, 1, rescale_mean_zero))
rownames(z_full) <- traits
png(paste0(dir,"b2_traits7_nocutoff_full.png"))
pheatmap(z_full)
dev.off()

breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b2_traits7_cutoff5_full.png"))
pheatmap(z_full,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()

# same thing on all of traits
h1 <- h1 %>% mutate(z=beta/se)
z <- tidyr::spread(h1 %>% select(trait,cluster,z),cluster,z)
traits <- z$trait
rownames(z) <- traits
z <- z %>% select(-trait)
png(paste0(dir,"b2_alltraits_nocutoff.png"))
pheatmap(z,show_rownames=F)
dev.off()
breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b2_alltraits_cutoff5.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,show_rownames=F)
dev.off()
## add in dropped cluster w/ rescale
z <- z %>% mutate(dropped=0)
# Function to rescale row values to have mean 0
rescale_mean_zero <- function(row) {
  mean_zero <- row - mean(row)
  return(mean_zero)
}
# Apply the rescaling function to each row of the dataframe
z_full <- t(apply(z, 1, rescale_mean_zero))
rownames(z_full) <- traits
png(paste0(dir,"b2_alltraits_nocutoff_full.png"))
pheatmap(z_full,show_rownames=F)
dev.off()

breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b2_alltraits_cutoff5_full.png"))
pheatmap(z_full,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,show_rownames=F)
dev.off()

####
## BLOCK 3
h1 <- data.table::fread("block3_9clusters.txt") %>% rename(trait=nimi)
filter(h1,p<0.005)

# how many traits have at least one cluster w/ p < 0.005
length(unique(filter(h1,p<0.005)$trait)) # 164 (of the 269 tested)

# filter to a small number of traits of interest (have >= 6 clusters w/ p<0.005)
head(as.data.frame(filter(h1,p<0.005) %>% group_by(trait) %>% count() %>% arrange(-n)),15)
trait3 <- (as.data.frame(filter(h1,p<0.005) %>% group_by(trait) %>% count() %>% arrange(-n)) %>% filter(n>=6))$trait
h <- filter(h1,trait %in% trait3)
h <- h %>% mutate(z=beta/se)

beta <- tidyr::spread(h %>% select(trait,cluster,beta),cluster,beta)

z <- tidyr::spread(h %>% select(trait,cluster,z),cluster,z)
traits <- z$trait
rownames(z) <- traits
z <- z %>% select(-trait)
png(paste0(dir,"b3_traits6_nocutoff.png"))
pheatmap(z)
dev.off()
breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b3_traits6_cutoff5.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()
breaksList = seq(-10, 10, by = 1)
png(paste0(dir,"b3_traits6_cutoff10.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()

## add in dropped cluster w/ rescale
z <- z %>% mutate(dropped=0)
# Function to rescale row values to have mean 0
rescale_mean_zero <- function(row) {
  mean_zero <- row - mean(row)
  return(mean_zero)
}
# Apply the rescaling function to each row of the dataframe
z_full <- t(apply(z, 1, rescale_mean_zero))
rownames(z_full) <- traits
png(paste0(dir,"b3_traits6_nocutoff_full.png"))
pheatmap(z_full)
dev.off()

breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b3_traits6_cutoff5_full.png"))
pheatmap(z_full,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)
dev.off()

# same thing on all of traits
h1 <- h1 %>% mutate(z=beta/se)
z <- tidyr::spread(h1 %>% select(trait,cluster,z),cluster,z)
traits <- z$trait
rownames(z) <- traits
z <- z %>% select(-trait)
png(paste0(dir,"b3_alltraits_nocutoff.png"))
pheatmap(z,show_rownames=F)
dev.off()
breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b3_alltraits_cutoff5.png"))
pheatmap(z,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,show_rownames=F)
dev.off()
## add in dropped cluster w/ rescale
z <- z %>% mutate(dropped=0)
# Function to rescale row values to have mean 0
rescale_mean_zero <- function(row) {
  mean_zero <- row - mean(row)
  return(mean_zero)
}
# Apply the rescaling function to each row of the dataframe
z_full <- t(apply(z, 1, rescale_mean_zero))
rownames(z_full) <- traits
png(paste0(dir,"b3_alltraits_nocutoff_full.png"))
pheatmap(z_full,show_rownames=F)
dev.off()

breaksList = seq(-5, 5, by = 1)
png(paste0(dir,"b3_alltraits_cutoff5_full.png"))
pheatmap(z_full,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,show_rownames=F)
dev.off()
