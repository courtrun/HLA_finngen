library(dplyr)

args = commandArgs(trailingOnly=TRUE)

b=commandArgs(TRUE)[1] # "/oak/stanford/groups/pritch/scratch.copy/groups/pritch/nasa/tools/ldsc/1000G_EUR_Phase3_plink/mymerge.bim"
s=commandArgs(TRUE)[2] # "/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/I9_CVD_HARD.gz"
o=commandArgs(TRUE)[3] # "/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hg19/I9_CVD_HARD.gz"

bim <- data.table::fread(b)
sumstats <- data.table::fread(s)

ss <- left_join(sumstats,bim %>% select(V1,V2,V4),by=c("rsids"="V2")) %>%
      rename(POS_hg38=pos,BP=V4,CHR_hg38=`#chrom`,CHR=V1,SNP=rsids,A1=alt,A2=ref)
ss <- ss %>% mutate(Z=beta/sebeta,N=412181)
ss <- filter(ss,!is.na(CHR))

write.table(ss,o, quote=F, sep="\t", row.names=F, col.names=T)

