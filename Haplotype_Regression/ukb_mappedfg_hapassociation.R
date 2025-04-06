
############ Identify overlapping snps between UKBB and FG
ml R/4.0
R
library(dplyr)

fg_full <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/allmaf01biallelicsnps.txt",header=T)

fg <- fg_full %>% filter(block=="1" & inhap=="yes")
uk <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/snp_list_maf01_hlablock1.txt.snplist",header=F)
overlapdf <- filter(uk,V1 %in% fg$rsids)
nrow(overlapdf)
data.table::fwrite(overlapdf, "/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb1_snps_overlap.txt", row.names = F, col.names=F, quote = F, sep = "\t")

fg <- fg_full %>% filter(block=="2" & inhap=="yes")
uk <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/snp_list_maf01_hlablock2.txt.snplist",header=F)
overlapdf <- filter(uk,V1 %in% fg$rsids)
nrow(overlapdf)
data.table::fwrite(overlapdf, "/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb2_snps_overlap.txt", row.names = F, col.names=F, quote = F, sep = "\t")

fg <- fg_full %>% filter(block=="3" & inhap=="yes")
uk <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/snp_list_maf01_hlablock3.txt.snplist",header=F)
overlapdf <- filter(uk,V1 %in% fg$rsids)
nrow(overlapdf)
data.table::fwrite(overlapdf, "/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb3_snps_overlap.txt", row.names = F, col.names=F, quote = F, sep = "\t")


############# finngen haplotypes for just the overlapped snps
# replace header with rsids then filter
library(dplyr)

an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/haplotypes/allmaf01biallelicsnps.txt") # file with annotated info

# block 1
a <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/supptable3_tab2_block1.csv") # full finngen hap
os <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb1_snps_overlap.txt",header=F) # overlap rsids to filter to
os_ids <- (filter(an,rsids %in% os$V1 & inhap=="yes") %>% arrange(pos))$V1 # filter to overlap
a_sub <- a %>% select(os_ids,total_doses,haplotype_group) # subset down to overlap
colnames(a_sub) <- c((filter(an,rsids %in% os$V1 & inhap=="yes") %>% arrange(pos))$rsids,"total_doses","haplotype_group") # update header name
data.table::fwrite(a_sub, "/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/fg_block1_overlap.txt", row.names = F, col.names=T, quote = F, sep = "\t")

a <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/supptable3_tab3_block2.csv")
os <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb2_snps_overlap.txt",header=F) # overlap rsids to filter to
os_ids <- (filter(an,rsids %in% os$V1 & inhap=="yes") %>% arrange(pos))$V1 # filter to overlap
a_sub <- a %>% select(os_ids,total_doses,haplotype_group) # subset down to overlap
colnames(a_sub) <- c((filter(an,rsids %in% os$V1 & inhap=="yes") %>% arrange(pos))$rsids,"total_doses","haplotype_group") # update header name
data.table::fwrite(a_sub, "/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/fg_block2_overlap.txt", row.names = F, col.names=T, quote = F, sep = "\t")

a <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/supptable3_tab4_block3.csv")
os <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb3_snps_overlap.txt",header=F) # overlap rsids to filter to
os_ids <- (filter(an,rsids %in% os$V1 & inhap=="yes") %>% arrange(pos))$V1 # filter to overlap
a_sub <- a %>% select(os_ids,total_doses,haplotype_group) # subset down to overlap
colnames(a_sub) <- c((filter(an,rsids %in% os$V1 & inhap=="yes") %>% arrange(pos))$rsids,"total_doses","haplotype_group") # update header name
data.table::fwrite(a_sub, "/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/fg_block3_overlap.txt", row.names = F, col.names=T, quote = F, sep = "\t")

########### EXTRACT, PROCESS UKBB HAPLOTYPES
#### Extract the genotypes as a phased vcf file at those n snps for just white uk biobank individuals, coded as 0/1
cd /oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data
ml biology
ml plink/2.0a2
plink2 --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --export vcf --extract hb1_snps_overlap.txt --output-missing-genotype 2 --out hb1_genotypes01_overlap

plink2 --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --export vcf --extract hb2_snps_overlap.txt --output-missing-genotype 2 --out hb2_genotypes01_overlap

plink2 --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --export vcf --extract hb3_snps_overlap.txt --output-missing-genotype 2 --out hb3_genotypes01_overlap

#### Process the vcf files to extract the phased data per person into haplotypes
############# MUNGE VCF, reformat to cleaned text file
ml R/4.0
R

library(dplyr)
blocks=1:3
nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"

for (bnum in blocks){
v <- data.table::fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb",bnum,"_genotypes01_overlap.vcf"))
data.table::fwrite(v, paste0(nsnp_dir,"hb",bnum,"_vcftodf_overlap.txt"), row.names = F, quote = F, sep = "\t")
}

############# SPLIT PHASING INTO SEPARATE COLUMNS
cd /oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/

### BLOCK 1
# Define files
input="hb1_vcftodf_overlap.txt"
output="hb1_vcftodf_phased_overlap.txt"

# Generate header
head -n 1 $input | tr '\t' '\n' | awk '{print (NR<10)?$0:$0"_a\t"$0"_b"}' | tr '\n' '\t' | sed 's/\t$/\n/' > $output

# Split phasing
awk 'BEGIN{ FS=OFS="\t" } {for (i=10; i<=NF; i++) {split($i, a, "|"); $i=a[1] OFS a[2]}} 1' $input  | tail -n +2 >> $output

### BLOCK 2
# Define files
input="hb2_vcftodf_overlap.txt"
output="hb2_vcftodf_phased_overlap.txt"

# Generate header
head -n 1 $input | tr '\t' '\n' | awk '{print (NR<10)?$0:$0"_a\t"$0"_b"}' | tr '\n' '\t' | sed 's/\t$/\n/' > $output

# Split phasing
awk 'BEGIN{ FS=OFS="\t" } {for (i=10; i<=NF; i++) {split($i, a, "|"); $i=a[1] OFS a[2]}} 1' $input  | tail -n +2 >> $output

### BLOCK 3
# Define files
input="hb3_vcftodf_overlap.txt"
output="hb3_vcftodf_phased_overlap.txt"

# Generate header
head -n 1 $input | tr '\t' '\n' | awk '{print (NR<10)?$0:$0"_a\t"$0"_b"}' | tr '\n' '\t' | sed 's/\t$/\n/' > $output

# Split phasing
awk 'BEGIN{ FS=OFS="\t" } {for (i=10; i<=NF; i++) {split($i, a, "|"); $i=a[1] OFS a[2]}} 1' $input  | tail -n +2 >> $output

############## assign ukb haps to nearest finngen hap groups then run hap association disease association
ml R/4.4.2
R
library(dplyr)
library(Rfast)
library(ggplot2)

nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"
blocks=c(1,2,3)

for (bnum in blocks){

  # load data
  v <- data.table::fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb",bnum,"_vcftodf_phased_overlap.txt")) # snps present in ukb and fg haps phased data

  # select "_a" phased genotypes
  selected_cols <- c(3,seq(10, ncol(v), by = 2)) # keep just the variant ID and all "_a" phased genotypes
  haplo1a <- v[, ..selected_cols]
  rownames(haplo1a) <- haplo1a$ID
  haplo1a <- haplo1a %>% select(-ID)
  snps <- rownames(haplo1a)
  haplo1at <- as.data.frame(t(haplo1a))
  colnames(haplo1at) <- snps

  # select "_b" phased genotypes
  selected_cols <- c(3,seq(11, ncol(v), by = 2)) # keep just the variant ID and all "_b" phased genotypes
  haplo1b <- v[, ..selected_cols]
  rownames(haplo1b) <- haplo1b$ID
  haplo1b <- haplo1b %>% select(-ID)
  snps <- rownames(haplo1b)
  haplo1bt <- as.data.frame(t(haplo1b))
  colnames(haplo1bt) <- snps

  # get IDs
  id_list <- gsub("_b","",rownames(haplo1bt))

  # combine phased genotypes
  haplo <- bind_rows(haplo1at,haplo1bt)
  print(dim(haplo)) # 2x number of people
  m <- distinct(haplo) # df with n snps as columns and 2x individuals as rows, keeping only distinct haplotypes
  print(dim(m)) # number unique haplotypes across all people's two copies
  mh <- tidyr::unite(m,hap,everything(),sep="") # create a new column "hap" that combines all other columns
  rownames(m) <- mh$hap # Add rownames with haplo info

  # Create dataframes with info in alt ways
  h <- tidyr::unite(haplo,hap,everything(),sep="") # create a new column "hap" that combines all other columns
  first <- data.frame(ID=paste0("ID ",c(1:(nrow(h)/2))),haplo1a=(h[c(1:(nrow(h)/2)),]),haplo1b=(h[c((nrow(h)/2+1):(nrow(h))),])) # haplotypes per person (one column per copy, one row per person)
  first$haplo1a <- as.character(first$haplo1a)
  first$haplo1b <- as.character(first$haplo1b)
  h$ID <- paste0("ID ",rep(c(1:(nrow(h)/2)),2))
  t1 <- h %>% group_by(hap) %>% count() %>% rename(total=n) # how many total doses each haplotype has (can be 2 from one individual)z
  t2 <- distinct(h) %>% group_by(hap) %>% count() %>% rename(unippl=n) # how many unique people have at least one dose of each haplotyep
  t3 <- filter(first,haplo1a==haplo1b) %>% group_by(haplo1a) %>% count() %>% rename(hap=haplo1a,twodose=n) # how many have 2 doses of each haplotype
  t4 <- data.frame(hap=c(filter(first,haplo1a!=haplo1b)$haplo1a,filter(first,haplo1a!=haplo1b)$haplo1b)) %>% group_by(hap) %>% count() %>% rename(onedose=n) # number one dose each
  t <- full_join(full_join(full_join(t1,t2,by=c("hap")),t3,by=c("hap")),t4,by=c("hap"))
  t[is.na(t)] <- 0
  mmh <- bind_cols(m,mh)
  m_full <- left_join(mmh,t,by=c("hap"))
  m_full$dendroname <- paste0("Hap: ",m_full$hap)

# filter to non-rare haplotypes
m_full <- m_full %>% filter(total>5) # >5 had 9879, >10 had 5914 remaining unique haplotypes

# cluster with finngen haps
fh <- data.table::fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/fg_block",bnum,"_overlap.txt"))
fh <- as.data.frame(fh)

# number of snps
nhapsnps <- ncol(fh) - 2 # subtract the total_doses, and haplotype_group columns

# Initialize
distances <- list()
min_distances <- list()
min_indices <- list()
jresults <- list()
all_dist_results <- list()

# loop through each row and calculate min euclidean distance
for (i in 1:nrow(m_full)) {
    distances <- Rfast::Dist(rbind(m_full[i, 1:nhapsnps], fh[, 1:nhapsnps]))[1, -1]
    min_distances[i] <- min(distances)
    min_indices[i] <- which.min(distances)
    print(i)
}

	min_dist_results <- data.frame(Row_in_ukb = 1:nrow(m_full),
	                      Min_Distance = unlist(min_distances),
	                      Closest_Row_in_fh = unlist(min_indices))
	#all_dist_results_comb <- bind_rows(all_dist_results)


# analyze
summary(min_dist_results$Min_Distance) # summary stats for the min dist assignments across ukb haplotypes, # min 0, median 2.45, mean 2.46, max 5.1
png(paste0("/home/users/courtrun/hist_min_dist_results_bnum",bnum,"_greater5dose.png"))
hist(min_dist_results$Min_Distance)
dev.off()

head(fh[(min_dist_results$Closest_Row_in_fh),ncol(fh)]) # finngen hap group assignments for ukb haps
data.frame(nearest_fgcluster=fh[(min_dist_results$Closest_Row_in_fh),ncol(fh)]) %>% group_by(nearest_fgcluster) %>% count() # nuber ukb indiv hap in each finngen hap group

# add assignments in
m_full$final_cluster <- fh[(min_dist_results$Closest_Row_in_fh),ncol(fh)]
m_full %>% group_by(final_cluster) %>% summarize(sum(total))

data.table::fwrite(m_full, paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t")

#######
### MUNGE FORMAT FOR REGRESSION
  #all <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01_overlap_greater5dose.txt")) # read in the above if not running script from start
  all <- m_full
  m_full <- 1
  l <- distinct(all %>% rename(Cluster=final_cluster,Haplotype=hap)) # long form of haplotypes, num individudals with that haplotype, cluster that haplotype belongs to
  temp<-distinct(left_join(first,l,by=c("haplo1a"="Haplotype")) %>% rename(cluster1a=Cluster))
  hap_an<-distinct(left_join(temp, l, by=c("haplo1b"="Haplotype")) %>% rename(cluster1b=Cluster))
  hap_an$UKB_ID <- id_list # add list of IDs from v
  hap_an <- filter(hap_an,!is.na(cluster1b)&!is.na(cluster1a)) # remove NA assignments if relevant
  nclust=max(hap_an$cluster1a,hap_an$cluster1b)
  for (n in 1:nclust){ # add one column per cluster and code 0, 1, or 2 for doses each person has of that cluster
    varname = paste0("cluster",n)
    hap_an[, varname] <- as.numeric(as.character(ifelse(hap_an$cluster1a==n&hap_an$cluster1b==n,2,ifelse(hap_an$cluster1a==n|hap_an$cluster1b==n,1,0))))
  }

  print(paste0("I've been here 567: block",bnum))

  # Make df of just ID and cluster belonging
  cassign <- hap_an %>% select(UKB_ID,matches("cluster")) %>% select(-cluster1a,-cluster1b)
  data.table::fwrite(cassign, paste0(nsnp_dir,"hb",bnum,"_clusterdoses_overlap_greater5dose_greater5dose.txt"), row.names = F, quote = F, sep = "\t")
  cassign <- cassign %>% select(-UKB_ID)

  # Calculate doses
  cdoses <- data.frame(
    Cluster = 1:ncol(cassign),
    cluster_onedose = colSums(apply(cassign, 2, function(x) x == 1)),
    cluster_twodose = colSums(apply(cassign, 2, function(x) x == 2)),
    cluster_unippl = colSums(apply(cassign, 2, function(x) x %in% c(1, 2)))
  ) %>% mutate(cluster_total=cluster_twodose*2+cluster_onedose)

  print(paste0("I've been here 678: block",bnum))
  data.table::fwrite(cdoses, paste0(nsnp_dir,"hb",bnum,"_clusterstats_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t")
  data.table::fwrite(l, paste0(nsnp_dir,"hb",bnum,"_haplostats_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t")

  # Add in covariate info
  cov_all <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe") %>%
    mutate(ID=paste0(`#FID`,"_",`IID`))
  mycols_cov <- c("ID","sex","age","PC1","PC2","PC3",
                  "PC4","PC5","PC6","PC7","PC7",
                  "PC8","PC9","PC10")
  cov_cov <- cov_all %>% select(matches(mycols_cov))

  # Add and recode disease data (currently 1,2,-9) to be 0/1 for no/yes has disease, set -9 to NA
  i <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe.info.tsv")
  mycols_traits <- unique(i$`#GBE_ID`)
  cov_traits <- cov_all %>% select(mycols_traits)

  # clear space
  cov_all <- 1

  # select just the disease coded columns and recode
  selected_columns <- cov_traits[, sapply(cov_traits, function(column) all(column %in% c(-9, 1, 2, NA)))] # select all columns with ONLY 1,2,-9 or NA as coding
  cov_traits <- cov_traits %>% select(colnames(cov_traits)[selected_columns])
  cov_traits[cov_traits==-9]<- NA # set missing
  cov_traits[cov_traits==1]<- 0 # set as no disease
  cov_traits[cov_traits==2]<- 1 # set as disease

  selected_traits <- unique(data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/traits_tested_ukbb.txt",fill=TRUE)$trait_ID)
  cov_traits <- cov_traits %>% select(selected_traits)

  data.table::fwrite(hap_an, paste0(nsnp_dir,"hb",bnum,"_temphapan_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t")

  # combine
  cov <- bind_cols(cov_cov,cov_traits) # combine covariate and disease info
  dat <- left_join(hap_an,cov, by=c("UKB_ID"="ID"))

  data.table::fwrite(dat, paste0(nsnp_dir,"hb",bnum,"_beforeregandleaveoutalldiseases_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t")

   #### Prepare to run regressions
  # Leave out the most frequent cluster
  leaveout = paste0("cluster",(cdoses %>% filter(cluster_total==max(cdoses$cluster_total)))$Cluster)
  print(leaveout) # cluster 5 for block 1
  keptclusters <- colnames(dat %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b,-leaveout))
  nclust <- length(keptclusters)

  data.table::fwrite(dat %>% select(-leaveout), paste0(nsnp_dir,"hb",bnum,"_preregressiondf_wleaveout_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t") # save df before running regression

  # remove clusters (and corresponding individuals) that have < 1000 total doses
	dat_full <- dat
	dat <- dat %>% filter(total.x > 1000 & total.y > 1000)
  nrow(dat %>% filter(total.x > 1000 & total.y > 1000)) # how many people remaining, block 1 254677, when filter to both > 5000 then leaves 206676
  oldnclust <- nclust
  keptclusters <- colnames(dat %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b,-leaveout))

# update list of kept clusters to remove those that have no assignments (commenting this out bc dont need since results are same since regress treats these already like they were removed
#zero_cols <- names(dat %>% select(matches("cluster")))[sapply(dat %>% select(matches("cluster")), function(col) all(col == 0))]
#keptclusters <- setdiff(keptclusters,zero_cols)

  # initialize data
  results_dt=data.table::data.table(trait=character(),beta=numeric(),se=numeric(),p=numeric())

  for (trait in colnames(cov_traits)){ # trait=mycols_traits[12]
    cov_traits %>% group_by(!!sym(trait)) %>% count()
    form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+sex+age+PC1+PC2+PC3+
      PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=dat)
    keptclusters <- (data.frame(coef=rownames(summary(regres)$coefficients)) %>% filter(grepl("cluster",coef)))$coef
    nclust <- length(keptclusters)
    print(summary(regres)$coefficients)
    beta<-summary(regres)$coefficients[2:(nclust+1),1] # select all cluster info (skip intercept and covariates)
    se<-summary(regres)$coefficients[2:(nclust+1),2]
    p<-summary(regres)$coefficients[2:(nclust+1),4]
    tmp_dt=data.table::data.table("trait"=trait,"beta"=beta,"se"=se,"p"=p)
    results_dt=data.table::rbindlist(list(results_dt,tmp_dt))
    print(trait)
  }

  # save regression results
  cluster<-gsub("cluster","",keptclusters)
  results_dt$cluster<-rep(cluster, length(colnames(cov_traits))) # results_dt$cluster<-cluster

  data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_afterregalldiseases_overlap_greater5dose.txt"), row.names = F, quote = F, sep = "\t")

  print(paste0("Done regression: block",bnum))

  }

  ############## Compare UKBB and FinnGen results
  library(dplyr)
library(ggplot2)

nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"

t <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/traits_tested_ukbb.txt",fill=TRUE)

## Load data and get stats
r1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_afterregalldiseases_overlap_greater5dose.txt")) %>%
	mutate(z=beta/se,block=1)
r2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_afterregalldiseases_overlap_greater5dose.txt")) %>%
	mutate(z=beta/se,block=2)
r3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_afterregalldiseases_overlap_greater5dose.txt")) %>%
	mutate(z=beta/se,block=3)
r <- bind_rows(r1,r2,r3)

traitmatch <- data.frame(
  uktrait = c( "HC151", "HC219", "BIN21068", "HC201", "HC1201",
              "HC430", "HC321", "BIN2443", "HC59", "HC136", "HC497",
              "HC78", "HC1211", "HC438", "HC1213"),
  fgtrait = c("G6_MS", "HYPOTHY_REIMB",
              "K11_COELIAC", "K11_IBD_STRICT", "L12_LUPUS",
              "RHEUMA_NOS", "M13_SJOGREN", "KELA_DIAB_INSUL_EXMORE",
              "THYROTOXICOSIS", "M13_SYSTCONNECT", "AB1_SEXUAL_TRANSMISSION",
              "M13_ARTHROPATHIES", "RHEUMA_SEROPOS_OTH", "K11_CHRONHEP",
				 "M13_PSORIARTH"),
  stringsAsFactors = FALSE
)

fgr <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/supptable4_tab1.csv")

# check if traits have at least one |Z|>4 in finngen
fgr %>% filter(traits %in% traitmatch$fgtrait) %>% filter(abs(Z_rescaled)>4) %>% group_by(traits) %>% count()

rfilt <- distinct(inner_join(r,traitmatch,by=c("trait"="uktrait")))
fgrfilt <- distinct(inner_join(fgr,traitmatch,by=c("traits"="fgtrait")))

b <- distinct(inner_join(rfilt,fgrfilt,by=c("trait"="uktrait","fgtrait"="traits","cluster"="hapgroup","block")))

z_cut=2
b  %>% filter(abs(Z_rescaled)>z_cut) %>% summarize(cor(z, Z_rescaled),n())
cor.test((b  %>% filter(abs(Z_rescaled)>z_cut))$z,(b  %>% filter(abs(Z_rescaled)>z_cut))$Z_rescaled)
b  %>% filter(abs(Z_rescaled)>z_cut)%>% group_by(trait) %>% summarize(cor(z, Z_rescaled),n())
z_cut=4
b  %>% filter(abs(Z_rescaled)>z_cut) %>% summarize(cor(z, Z_rescaled),n())
cor.test((b  %>% filter(abs(Z_rescaled)>z_cut))$z,(b  %>% filter(abs(Z_rescaled)>z_cut))$Z_rescaled)
z_cut=0
b  %>% filter(abs(Z_rescaled)>z_cut) %>% summarize(cor(z, Z_rescaled),n())
cor.test((b  %>% filter(abs(Z_rescaled)>z_cut))$z,(b  %>% filter(abs(Z_rescaled)>z_cut))$Z_rescaled)

nrow(b  %>% filter(abs(Z_rescaled)>4)) #  associations that have |Z|>4 in FinnGenn
nrow(b  %>% filter(abs(Z_rescaled)>4 & sign(z)==sign(Z_rescaled))) # same sign
nrow(b  %>% filter(abs(Z_rescaled)>4 & abs(z)>1.645)) # |Z| > 1.645 in UKB and |Z|>4 in finngen
nrow(b  %>% filter(abs(Z_rescaled)>4 & abs(z)>1.645 & sign(z)==sign(Z_rescaled))) # |Z| > 1.645 in UKB and |Z|>4 in finngen and same sign

nrow(traitmatch) # how many traits included

# view traits
left_join(b,t,by=c("trait"="trait_ID")) %>% select(trait,fgtrait,GBE_NAME) %>% distinct()

# checking individual association replications
b  %>% filter(abs(Z_rescaled)>4 & abs(z)>2)

# making table to include in zenodo or supplement
bclean <- b %>% select(block,cluster,fgtrait,Z_rescaled,trait,z) %>%
	 rename(hapgroup=cluster,finngen_trait=fgtrait,
					 finngen_Z=Z_rescaled,ukb_trait=trait,ukb_Z=z) %>% arrange(block,finngen_trait)
data.table::fwrite(bclean, paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/ukbreplicationresults_comparefinngen.txt"), row.names = F, quote = F, sep = "\t")

binom.test(10,33,0.05)

ggplot(bclean,aes(finngen_Z,ukb_Z))+
    geom_point()+
    geom_smooth(method = "lm", se = TRUE, color = "blue") + theme_classic()+
    labs(x="FinnGen Haplotype Association Z-scores",
        y="UKBB Haplotype Association Z-scores")

########## Pearson's correlations
library(dplyr)
library(ggplot2)

nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"

t <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/traits_tested_ukbb.txt",fill=TRUE)

## Load data and get stats
r1 <- data.table::fread(paste0(nsnp_dir,"hb",1,"_afterregalldiseases_overlap_greater5dose.txt")) %>%
	mutate(z=beta/se,block=1)
r2 <- data.table::fread(paste0(nsnp_dir,"hb",2,"_afterregalldiseases_overlap_greater5dose.txt")) %>%
	mutate(z=beta/se,block=2)
r3 <- data.table::fread(paste0(nsnp_dir,"hb",3,"_afterregalldiseases_overlap_greater5dose.txt")) %>%
	mutate(z=beta/se,block=3)
r <- bind_rows(r1,r2,r3)

traitmatch <- data.frame(
  uktrait = c( "HC151", "HC219", "BIN21068", "HC201", "HC1201",
              "HC430", "HC321", "BIN2443", "HC59", "HC136", "HC497",
              "HC78", "HC1211", "HC438", "HC1213"),
  fgtrait = c("G6_MS", "HYPOTHY_REIMB",
              "K11_COELIAC", "K11_IBD_STRICT", "L12_LUPUS",
              "RHEUMA_NOS", "M13_SJOGREN", "KELA_DIAB_INSUL_EXMORE",
              "THYROTOXICOSIS", "M13_SYSTCONNECT", "AB1_SEXUAL_TRANSMISSION",
              "M13_ARTHROPATHIES", "RHEUMA_SEROPOS_OTH", "K11_CHRONHEP",
             "M13_PSORIARTH"),
  stringsAsFactors = FALSE
)

fgr <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/supptable4_tab1.csv")

rfilt <- distinct(inner_join(r,traitmatch,by=c("trait"="uktrait")))
fgrfilt <- distinct(inner_join(fgr,traitmatch,by=c("traits"="fgtrait")))

b <- distinct(inner_join(rfilt,fgrfilt,by=c("trait"="uktrait","fgtrait"="traits","cluster"="hapgroup","block")))

trait_list <- unique(b$trait)

# For loop to iterate through each unique pairwise combination
count=0
cor_results <- list()
for (ibnum in 1:3) { # loop through each block
for (i in 1:(length(trait_list) - 1)) {
  for (j in (i + 1):length(trait_list)) {
	  count=count+1
    pair <- paste0(trait_list[i],"-", trait_list[j])
		t1z <- filter(b,block==ibnum & trait==trait_list[i])$z
		t2z <- filter(b,block==ibnum & trait==trait_list[j])$z
		cor.test(t1z,t2z)
		cor_results[[count]] <- data.frame(t1=trait_list[i],t2=trait_list[j],pair=pair,
					cor=cor.test(t1z,t2z)$estimate,p=cor.test(t1z,t2z)$p.value,block=ibnum)
    print(pair)
  }}}
  cor_comb <- bind_rows(cor_results)
cor_comb <- left_join(left_join(cor_comb,traitmatch,by=c("t1"="uktrait")) %>%
 rename(fgtrait1=fgtrait),traitmatch,by=c("t2"="uktrait")) %>% rename(fgtrait2=fgtrait) %>% # annotate
 mutate(fgtp1=paste0(fgtrait1,"-",fgtrait2),fgtp2=paste0(fgtrait2,"-",fgtrait1))

# load finngen data
fgcor <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/traitpair_correlation_ukbtraits.txt") # from finngen
fgcor <- fgcor %>% filter(traitpair %in% unique(c(cor_comb$fgtp1,cor_comb$fgtp2))) # filter to only the relevant trait pairs; has 5n rows bc 5 measures
gctraits <- unique((fgcor %>% filter(measure=="gc"&cor>0.3))$traitpair)
gctraitsallpos <- (fgcor %>% filter(measure %in% c("hapregb1","hapregb2","hapregb3") &traitpair %in% gctraits) %>%
			group_by(traitpair) %>% summarize(allpos=all(cor>0.3)) %>% filter(allpos==TRUE))$traitpair
gctraitsanyneg <-  (fgcor %>% filter(measure %in% c("hapregb1","hapregb2","hapregb3") &traitpair %in% gctraits) %>%
			group_by(traitpair)  %>% summarize(anyneg=any(cor < (-0.3))) %>% filter(anyneg==TRUE))$traitpair

gc3_uk <- filter(cor_comb,fgtp1 %in% gctraits | fgtp2 %in% gctraits) # has 3n rows bc 3 blocks
nrow(gc3_uk)/3 # divide by 3 because 3 blocks

gc3_uk_allpos <- (gc3_uk %>% group_by(pair) %>% summarize(allpos=all(cor>0.3)) %>% filter(allpos==TRUE))$pair
gc3_uk_allpos_fgids <- unique(c((cor_comb %>% filter(pair %in% gc3_uk_allpos))$fgtp1,(cor_comb %>% filter(pair %in% gc3_uk_allpos))$fgtp2)) # note this is 2x bc putting in both orders depending on how id match
cor_comb %>% filter(pair %in% gc3_uk_allpos & (fgtp1 %in% gctraitsallpos | fgtp2 %in% gctraitsallpos))
intersect(gctraitsallpos,gc3_uk_allpos_fgids) # trait pairs that are allpos in fg and ukb
binom.test(5,12,0.05)

gc3_uk_anyneg <- (gc3_uk %>% group_by(pair) %>% summarize(anyneg=any(cor < (-0.3))) %>% filter(anyneg==TRUE))$pair
gc3_uk_anyneg_fgids <- unique(c((cor_comb %>% filter(pair %in% gc3_uk_anyneg))$fgtp1,(cor_comb %>% filter(pair %in% gc3_uk_anyneg))$fgtp2))
intersect(gctraitsanyneg,gc3_uk_anyneg_fgids) # trait pairs that are anyneg in fg and ukb
binom.test(1,5,0.05)

# match together and write table, could add to supplement or zenodo
cor_comb2 <- cor_comb
cor_comb2 <- cor_comb2 %>% rename(fgtraitpair=fgtp2) %>% select(-fgtp1)
cor_comb1 <- cor_comb %>% rename(fgtraitpair=fgtp1) %>% select(-fgtp2)
fgcor1 <- left_join(fgcor %>% filter(measure %in% c("hapregb1","hapregb2","hapregb3")),
					fgcor %>% filter(measure=="gc"),by=c("traitpair")) %>%
			rename(fg_cor=cor.x,fg_p=p.x,gc_cor=cor.y,gc_p=p.y,block=measure.x) %>% select(-selow.x,-sehigh.x,-selow.y,-sehigh.y,-measure.y)
fgcor1$block <- as.numeric(gsub("hapregb","",fgcor1$block))
cor_comb_all <- inner_join(bind_rows(cor_comb1,cor_comb2),fgcor1,by=c("fgtraitpair"="traitpair","block")) %>%
rename(ukbtrait1=t1,ukbtrait2=t2,ukbtraitpair=pair,ukb_cor=cor,ukb_p=p)

data.table::fwrite(cor_comb_all, paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/extraoverlapfiles/ukb_comparefinngen_pearsonscor.txt"), row.names = F, quote = F, sep = "\t")
