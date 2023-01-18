# Courtney Smith and Satu Strausz - HLA project - Make Matrix
# Script started 1-18-2023
# Goal of script: Make matrix of clumped hits by traits
# Called in snakemake script by rule make_matrix

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

# Combine hits for all traits
my.list <- list()
for (i in 1:(length(args)-2)){
Trait_name=gsub("_filt428hits.gz","",gsub(".*filtered/","",commandArgs(TRUE)[i]))
temp <- data.table::fread(commandArgs(TRUE)[i],fill=TRUE)
colnames(temp) <- c("chr","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
my.list[[i]] <- temp %>% mutate(TraitName=Trait_name,ID=paste0("chr",chr,"_",pos,"_",ref,"_",alt)) %>% select(ID,rsids,pval,beta,sebeta,TraitName) %>%
 mutate(pval=as.numeric(pval))
}
trait_long <- bind_rows(my.list)

library(data.table)

trait_long$pval <- as.numeric(trait_long$pval)
trait_long$rsids <- ifelse(trait_long$rsids=="",trait_long$ID,trait_long$rsids)
write.table(trait_long,commandArgs(TRUE)[length(args)-1], quote=F, sep="\t", row.names=F, col.names=T)

# Long to wide form
trait_wide <- dcast(as.data.table(trait_long), ID + rsids ~ paste0("_", TraitName), fun=mean,value.var = c("beta","sebeta","pval"),sep='')
write.table(trait_wide,commandArgs(TRUE)[length(args)], quote=F, sep="\t", row.names=F, col.names=T)
