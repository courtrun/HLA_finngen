# Courtney Smith
# Created Feb 2024
# Goal of script: compare the traits found associated by SNP vs haplotype vs allele methods

library(dplyr)

hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=paste0(hla_dir,"snps",1000,"/")

# Set sig thresholds
sigthres1 = 6.7e-6 # for hap reg and per block allele reg; 0.05/(2459 traits *3 blocks)
sigthres2 = 2e-7 # for indiv allele reg; 0.05/(2459 traits * 98 alleles tested separately)

# Total number of significant associations for the haplotype regression analysis with all 2459 traits, not rescaled
results_hap1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres1) # from hapclusterandregress_alltraits.R
h1_traits <- unique(results_hap1$nimi)
results_hap2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres1) # from hapclusterandregress_alltraits.R
h2_traits <- unique(results_hap2$nimi)
results_hap3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres1) # from hapclusterandregress_alltraits.R
h3_traits <- unique(results_hap3$nimi)
h_traits <- unique(c(h1_traits,h2_traits,h3_traits))
results_haph <- bind_rows(results_hap1 %>% mutate(block=1),results_hap2 %>% mutate(block=2),results_hap3 %>% mutate(block=3))
nrow(results_haph %>% filter(p<6.7e-6)) # 1750
length(unique((results_haph %>% filter(p<6.7e-6))$nimi)) # 209
nrow(distinct((results_haph %>% filter(p<6.7e-6)) %>% select(cluster,block))) # 70
results_haph_all_traits <- results_haph

# allele results for one at a time across whole region all traits
results_allele <- data.table::fread(paste0(nsnp_dir,"alleles_indiv_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres2) # from allele_regress.R
a_traits <- unique(results_allele$nimi)
nrow(results_allele %>% filter(p<2e-7)) # 4835
length(unique((results_allele %>% filter(p<2e-7))$nimi)) # 377
nrow(distinct((results_allele %>% filter(p<2e-7)) %>% select(allele))) # 94

# Total number of significant associations for the haplotype regression analysis with the 272 HLA-associated traits as ascertained by SNP analysis
melted_df1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_heatmapplotted.txt"))
melted_df2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_heatmapplotted.txt"))
melted_df3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_heatmapplotted.txt"))
results_haph <- bind_rows(melted_df1 %>% mutate(block=1),melted_df2 %>% mutate(block=2),melted_df3 %>% mutate(block=3))
nrow(results_haph %>% filter(abs(value)>4))
nrow(results_haph %>% filter(abs(value)>4) %>% select(traits) %>% distinct()) # 64
nrow(results_haph %>% filter(abs(value)>4) %>% select(variable,block) %>% distinct()) # 65

# allele results for all within a given block together for all traits
results_allele_combA <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_blockA_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres1) # from allele_regress.R
acombA_traits <- unique(results_allele_combA$nimi)
results_allele_combC <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_blockC_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres1) # from allele_regress.R
acombC_traits <- unique(results_allele_combC$nimi)
results_allele_combD <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_blockD_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% filter(p < sigthres1) # from allele_regress.R
acombD_traits <- unique(results_allele_combD$nimi)
results_allele_allblocks_comb <- bind_rows(results_allele_combA %>% mutate(block=1),results_allele_combC %>% mutate(block=2),results_allele_combD %>% mutate(block=3))
a_block <- unique(results_allele_allblocks_comb$nimi)
nrow(results_allele_allblocks_comb %>% filter(p<6.7e-6)) # 2526
nrow(results_allele_allblocks_comb %>% filter(p<6.7e-6) %>% select(nimi) %>% distinct()) # 323
nrow(results_allele_allblocks_comb %>% filter(p<6.7e-6) %>% select(allele,block) %>% distinct()) # 60

# snp association results
results_snp_uncor <- data.table::fread("/finngen/green/courtrun/R10_condit_results.txt")
results_snp_cor <- data.table::fread("/finngen/green/courtrun/R10_correlated_traits_condtionals.txt")
snp_trait <- unique(c(results_snp_uncor$trait,results_snp_cor$trait)) # already filtered to the sig associated traits in at least one snp reg 1e-6
results_snp <- bind_rows(results_snp_cor,results_snp_uncor)
nrow(results_snp %>% filter(pval < 1e-6)) # 1064, sig associations in all traits
nrow(results_snp %>% filter(pval < 1e-6) %>% select(trait) %>% distinct) # 572, sig traits when considering all traits
nrow(results_snp %>% filter(pval < 1e-6)%>% select(ID) %>% distinct) # 687, sig variant associations when considering all traits
length(unique(results_snp_uncor$ID)) # 428 unique variants sig in the 272 non-redundant traits
nrow(results_snp_uncor) # 540 sig associations in the 272 non-redundant traits
length(unique(results_snp_uncor$trait)) # 272 unique traits sig in the 272 non-redundant traits

# Haplotype regressions adjusted for the alleles, performed on all traits
comb1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults_adj_all_alleles_inblock_v3.txt"))
comb2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults_adj_all_alleles_inblock_v3.txt"))
comb3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults_adj_all_alleles_inblock_v3.txt"))
comb <- bind_rows(comb1 %>% mutate(block=1),comb2 %>% mutate(block=2),comb3 %>% mutate(block=3)) %>% mutate(z_adj=beta/se)
nrow(filter(comb,abs(z_adj)>4&grepl("cluster",cluster))) # num sig associations in adj regression
length(unique((filter(comb,abs(z_adj)>4&grepl("cluster",cluster)))$nimi)) # num sig traits in adj regression
nrow((filter(comb,abs(z_adj)>4&grepl("cluster",cluster))) %>% select(cluster,block) %>% distinct()) # num sig hap groups in adj regression

###### ANALYSIS
# Total num sig traits in each
length(unique(h_traits))
length(unique(a_traits))
length(unique(a_block))
length(snp_trait)

# Total num sig traits in hap but not other methods
length(setdiff(h_traits,a_traits))
length(setdiff(h_traits,a_block))
length(setdiff(h_traits,snp_trait))
length(setdiff(h_traits,unique(c(a_traits,snp_trait))))
length(setdiff(h_traits,unique(c(a_block,snp_trait))))
length(setdiff(h_traits,unique(c(a_traits,a_block,snp_trait))))

# Dive into examples of the above
#filter(results_hap,nimi=="D3_PURPURA_AND3_OTHER_HAEMORRHAGIC")
#temp <- filter(results_hap,nimi %in% (setdiff(h_traits,unique(c(a_traits,a_block,snp_trait)))))

# Qualitative comparison of top associations for allele indv vs allele per blozck
temp1 <- head(results_allele_allblocks_comb %>% arrange(-abs(Z)),500)$nimi
temp2 <-head(results_allele %>% arrange(-abs(Z)),500)$nimi
length(unique(temp1))
length(unique(temp2))
length(unique(c(unique(temp1),unique(temp2))))
length(intersect(unique(temp1),unique(temp2)))

#### All together
# Number associations across methods
nrow(results_haph_all_traits)+nrow(results_allele)+nrow(results_snp) #only did approach 1 allele in this bc otherwise redundant #+results_allele_allblocks_comb
length(unique(c(results_haph_all_traits$nimi,results_allele$nimi,results_snp$trait,results_allele_allblocks_comb$nimi))) # included both allele methods bc not additive across double counting same traits

#### Saving SUPPLEMENTAL TABLE FILES
# Haplotype group regression across blocks
results_hap1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se)
results_hap2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) 
results_hap3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) 
results_haph <- bind_rows(results_hap1 %>% mutate(block=1),results_hap2 %>% mutate(block=2),results_hap3 %>% mutate(block=3)) %>%
  rename(trait=nimi,haplotype_group=cluster)
data.table::fwrite(results_haph, paste0(nsnp_dir,"all_hapgroup_reg_results.txt"), row.names = F, quote = F, sep = "\t")
results_haph_sig <- results_haph %>% filter(p<6.7e-6)
data.table::fwrite(results_haph_sig, paste0(nsnp_dir,"all_hapgroup_reg_results_sig.txt"), row.names = F, quote = F, sep = "\t")

# SNP Conditional analysis results
results_snp_uncor <- data.table::fread("/finngen/green/courtrun/R10_condit_results.txt")
results_snp_cor <- data.table::fread("/finngen/green/courtrun/R10_correlated_traits_condtionals.txt")
results_snp <- bind_rows(results_snp_cor,results_snp_uncor) %>% mutate(Z=beta/se)
data.table::fwrite(results_snp, paste0(nsnp_dir,"condit_snp_reg_results_sig.txt"), row.names = F, quote = F, sep = "\t") # just has sig associations

### Allele regression results
# allele results for one at a time across whole region all traits
results_allele <- data.table::fread(paste0(nsnp_dir,"alleles_indiv_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) %>% rename(trait=nimi)
data.table::fwrite(results_allele, paste0(nsnp_dir,"alleles_indiv_regresults_alltraits_supp.txt"), row.names = F, quote = F, sep = "\t")
results_allele_sig <- results_allele %>% filter(p<2e-7)
data.table::fwrite(results_allele_sig, paste0(nsnp_dir,"alleles_indiv_regresults_alltraits_sig_supp.txt"), row.names = F, quote = F, sep = "\t")

# allele results for all within a given block together for all traits
results_allele_combA <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_blockA_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se)
results_allele_combC <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_blockC_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se)
results_allele_combD <- data.table::fread(paste0(nsnp_dir,"alleles_vifindep_blockD_regresults_alltraits.txt")) %>% 
  mutate(Z=beta/se) 
results_allele_allblocks_comb <- bind_rows(results_allele_combA %>% mutate(block=1),results_allele_combC %>% mutate(block=2),results_allele_combD %>% mutate(block=3)) %>%
  rename(trait=nimi)
data.table::fwrite(results_allele_allblocks_comb, paste0(nsnp_dir,"alleles_vifindep_blocks_regresults_alltraits_supp.txt"), row.names = F, quote = F, sep = "\t")
results_allele_allblocks_comb_sig <- results_allele_allblocks_comb %>% filter(p<6.7e-6)
data.table::fwrite(results_allele_allblocks_comb_sig, paste0(nsnp_dir,"alleles_vifindep_blocks_regresults_alltraits_sig_supp.txt"), row.names = F, quote = F, sep = "\t")
