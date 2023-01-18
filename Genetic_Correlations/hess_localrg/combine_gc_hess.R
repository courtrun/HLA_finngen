# Created by Courtney Smith - HLA - LDSC Gen Cor / HESS Local Gen Cor
# Started on 1-11-2023
# Goal of script: Create file with combined info from LDSC genetic correlation pairs for hess step 3
# Called by snakefile in step prep_for_hess3b

library(dplyr)

c <- data.table::fread(commandArgs(TRUE)[1]) # /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/combined_gen_cor_sumstats.txt

c$p1 <- gsub(".*munged/","",c$p1)
c$p1 <- gsub(".sumstats.*","",c$p1)
c$p2 <- gsub(".*munged/","",c$p2)
c$p2 <- gsub(".sumstats.*","",c$p2)

c <- c %>% select(p1,p2,gcov_int) %>% rename(pheno_cor=gcov_int) # note: can just use gcov_int bc sample size are same and full overlap
c <- c %>% mutate(pheno_cor_bound = ifelse(pheno_cor > 1,1,ifelse(pheno_cor < -1,1,pheno_cor)))

OUTPUT_FILE=commandArgs(TRUE)[2]
write.table(c, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)
