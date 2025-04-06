### Courtney Smith - Enrichment Analysis - Prep genome sumstats hits files
### Started on 12-7-2023
### Goal: Plot number hits in genome and HLA region across all traits in bp blocks

library(dplyr)
library(ggplot2)
library(RColorBrewer)

a <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_gwas_hits/filt/all_hits_sumstats.txt",header=F)
colnames(a) <- c("trait","ID","CHR","BP","ref","alt","BETA","SE","P","AF_alt")

# filter to just pruned trait list
ptl <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/prunedalltraitsforenrichment.txt")
a <- filter(a,trait %in% ptl$traits)

# Define the size of the buckets (1 million in this case)
bucket_size <- 100000 # 100kb, 50kb, 500kb, 200kb
num_in_kb <- bucket_size/1000

# Create a new column representing the bucket for "BP"
a <- a %>% arrange(CHR,BP) %>%
  mutate(BP_bucket = floor(BP / bucket_size) * bucket_size)

# Group by "chrom" and "BP_bucket" and count the number of rows in each group
result <- a %>%
  group_by(CHR, BP_bucket) %>%
  summarise(count = n())

########## Add in HLA
# correlated traits
c <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_correlated_traits_condtionals.txt")
c <- distinct(c %>% select(trait,rsids,pos))
# uncorrelated traits
u <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_condit_results.txt")
u <- distinct(u %>% select(trait,rsids,pos))
h <- bind_rows(c,u) %>% arrange(pos)

# filter to just pruned traits
h <- filter(h,trait %in% ptl$traits)

nhla <- nrow(h)
h <- h %>% arrange(pos) %>%
  mutate(BP_bucket = floor(pos / bucket_size) * bucket_size)
# Group by "chrom" and "BP_bucket" and count the number of rows in each group
hresult <- h %>%
  group_by(BP_bucket) %>%
  summarise(count = n()) %>% mutate(CHR=6) %>% select(CHR,BP_bucket,count)
result <- bind_rows(result,hresult) %>% arrange(CHR,BP_bucket)

# Add continuous BP
#acont <- result %>% mutate(contBP=cumsum(BP_bucket))
result$order <- c(1:nrow(result))
result$CHR_color <- ifelse((result$CHR %% 2) == 0,"even","odd")

result$CHR <- as.character(result$CHR)
result <- bind_rows(filter(result,CHR!="23"),filter(result,CHR=="23") %>% mutate(CHR="X"))

# Set chromosome labeling
axis.set <- result %>%
  group_by(CHR) %>%
  summarize(center = (max(order) + min(order)) / 2)

# Plotting
ggplot(result,aes(x=order,y=count,color=as.factor(CHR_color)))+
  geom_point(size=5)+
   scale_color_brewer(palette="Dark2")+ guides(color = "none")+
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) + # label the chromosome number in the center of each chromosome block
  theme_classic()+labs(x="Chromosome",y=paste0("Number of Hits Across\nAll Traits in ",num_in_kb,"kb blocks")) +
  theme(text=element_text(size = 28),axis.text.x = element_text(size = 24))
ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/projects/hla/figures/enrichment_manhattan_hitsacrosstraits",bucket_size,"_prunedtraits.png"),width=30,height=8)

print(result %>% arrange(-count) %>% head)
summary(result$count) # sumstats for avg/min/median/max total finemapped hits across all traits in each bin across all bins
