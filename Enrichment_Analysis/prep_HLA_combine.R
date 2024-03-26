### Courtney Smith - Enrichment Analysis - Prep HLA and combine hits files
### Started on 3-5-2023
### Goal: Make a file that has enrichment in HLA for traits

library(dplyr)

# correlated traits
c <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_correlated_traits_condtionals.txt")
c <- distinct(c %>% select(trait,rsids))

# uncorrelated traits
u <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_condit_results.txt")
u <- distinct(u %>% select(trait,rsids))

h <- bind_rows(c,u) %>% arrange(trait)
hits <- h %>% group_by(trait) %>% count()

# compare to genome hits
g <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_gwas_hits/filt/all_hits.txt")
hits <- full_join(hits %>% rename(HLA_hits=n),g,by=c("trait"))
hits[is.na(hits)] <- 0

HLA_SNP=41234
genome_SNP=9727032

hits <-  hits %>% mutate(HLA_SNP=41234,genome_SNP=9727032)

# Across all traits
(mean(hits$HLA_hits)/HLA_SNP)/(mean(hits$non_HLA_hits)/genome_SNP) # 18.0
mean(hits$HLA_hits) # 1.03
mean(hits$non_HLA_hits) # 13.5

### Add in finngen annotations
an <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/data/FINNGEN_ENDPOINTS_DF10_Final_2022-05-16_public.txt")
an_long <- tidyr::separate_rows(an,"TAGS",sep=",") %>% select(TAGS,NAME)# split into multiple rows
an <- left_join(hits,an_long,by=c("trait"="NAME"))

# Add labels that the trait starts with for those that didn't automatically get labeled
unlab <- filter(an,is.na(TAGS))
unlab_split <- unlab %>% tidyr::separate(trait,sep="_",into=c("group","other")) %>% select(group)
unlab$TAGS <- unlab_split$group

t <- bind_rows(filter(an,!is.na(TAGS)),unlab)
t$TAGS <- gsub("#","",t$TAGS)
t <- t %>% rename(group=TAGS)

# then manually annotated and grouped, loading file below
tg <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/Finngen_group_trait_annotations.tsv") %>%
mutate(plotnames=ifelse(combine=="",group_name,supergroup_name))

t <- left_join(t,tg,by=c("group")) %>% rename(fg_group=group,group=plotnames)
OUTPUT_FILE="/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrich_trait_group_annotations.txt"
write.table(t, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

t <- distinct(t %>% select(-fg_group,-group_name)) # want to only count each trait one time per group so remove within group duplicates

tm <- t %>% group_by(group) %>% summarize(sd_HLA_hits=sd(HLA_hits),sd_non_HLA_hits=sd(non_HLA_hits),HLA_hits=mean(HLA_hits),non_HLA_hits=mean(non_HLA_hits))
tm <- left_join(tm,t %>% group_by(group) %>% count(),by=c("group"))
tm <-  tm %>% mutate(HLA_SNP=41234,genome_SNP=9727032) %>% mutate(enrich=(HLA_hits/HLA_SNP)/(non_HLA_hits/genome_SNP))

# Were any not enriched in HLA?
tm %>% filter(enrich<1)
e <- tm %>% select(group,enrich,n) %>% filter(n>5) # %>% filter(n>=20)
e <- left_join(e,tm,by=c("group","n","enrich"))

OUTPUT_FILE="/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrichHLA_bygroup_finngenan.txt"
write.table(e, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)
