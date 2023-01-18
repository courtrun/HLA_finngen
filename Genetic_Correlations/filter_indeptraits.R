# Courtney Smith - BOLT-LMM - Snakemake
# Script started 12-14-2022
# Goal of script: Filter freeze 10 traits to have no traits that are in rg > 0.95 with any other

library(dplyr)

r <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/combined_gen_cor_sumstats.txt")
r <- r %>% mutate(abs=abs(rg))
r <- na.omit(r %>% filter(!p1==p2))
r$abs <- ifelse(r$abs<0.01,0.01,r$abs) # prevent R from reading small values as characters
r$abs <- as.numeric(r$abs)

df_wide <- r %>% select(p1,p2,abs) %>%
  tidyr::pivot_wider(names_from = p2, values_from = abs)

# First identify traits that have abs(rg) < 0.95 with all traits
df_wide2 <- df_wide
df_wide2$max <- apply(df_wide[,-1],1,max,na.rm=TRUE)
df_wide2 <- df_wide2 %>% select(p1,max)
keep1 <- unique(filter(df_wide2,max<0.95)$p1)

# For the remaining traits, add trait one by one to list and other add traits not correlated with those
df_wide <- df_wide %>% filter(!(p1 %in% keep1))

keep2 <- list() # list of kept traits
removedtraits <- list() # list of traits that are correlated with a kept trait

for (item in df_wide$p1) {
if (item %in% removedtraits){
next}
correlated <- filter(r,(p1==item | p2==item)&abs>0.95) # identify all traits correlated with this trait
removetraits <- setdiff(unique(c(correlated$p1,correlated$p2)),item) # remove others
removedtraits <- c(removedtraits, removetraits) # add to list of removed traits
keep2 <- c(keep2, item) # add to list of kept traits
}

# Combine
keep <- unique(c(keep1,keep2))

# Convert to dataframe and add names
keep <- data.frame(V1=unlist(keep))
keep$V2 <- gsub(".sumstats.gz","",gsub(".*./","",as.character(keep$V1)))

# Write output
OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/independent_traits.tsv"
write.table(keep, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

### CHECKING to make sure it worked
# Remake original df_wide
df_wide <- r %>% select(p1,p2,abs) %>%
  tidyr::pivot_wider(names_from = p2, values_from = abs)

# For traits not kept, see if have a trait correlated w/ it in kept list
setdiff(df_wide$p1,temp$p1)
filter(temp,p1 %in% rownames(filter(data.frame(V1=t(filter(df_wide,p1=="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/munged/SLE_OTH.sumstats.gz"))),V1>0.95)))

# For traits kept, make sure none of its correlated were also kept
temp$p1
filter(temp,p1 %in% rownames(filter(data.frame(V1=t(filter(df_wide,p1=="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/munged/RX_CROHN_2NDLINE.sumstats.gz"))),V1>0.95)))
