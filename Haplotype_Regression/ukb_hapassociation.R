######### Get haplotype data for each block
# Load tools and environment
ml biology
ml plink/2.0a2
cd /oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data

##### Extract SNPs with MAF > 1% within HLA region, using the pfiles from the phased haplotypes (non-imputed data) from ukb24983 within each block
## Set parameters
MAF=0.01
OUT_FILE1="snp_list_maf01_hlablock1.txt"
OUT_FILE2="snp_list_maf01_hlablock2.txt"
OUT_FILE3="snp_list_maf01_hlablock3.txt"

## BLOCK 1
plink2 \
  --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --chr 6 \
  --from-bp 29590597 --to-bp 30013393 --maf $MAF \
  --write-snplist \
  --out $OUT_FILE1

## BLOCK 2
plink2 \
  --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --chr 6 \
  --from-bp 31136575 --to-bp 31578848 --maf $MAF \
  --write-snplist \
  --out $OUT_FILE2

## BLOCK 3
plink2 \
  --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --chr 6 \
  --from-bp 32062687 --to-bp 32814902 --maf $MAF \
  --write-snplist \
  --out $OUT_FILE3

##### Randomly sample to n snps in each block
ml R/4.0
R

library(dplyr)
library(data.table)
nsnp=1000
nsnp2=500
blocks=1:3

for (bnum in blocks){
  # Load in snps in the block, randomly subset
  #hb_all<-fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/snp_list_maf01_hlablock",bnum,".txt.snplist"), header = F) # all snps in the block w/ MAF > 1%
  hb_all_snp<-fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/snp_list_maf01_hlablock",bnum,".txt.snplist"), header = F) # all snps in the block w/ MAF > 1%
  #hb_all_snp <- hb_all[grep("^6:\\d+_[A-Z]_[A-Z]$", hb_all$V1), , drop = FALSE] # make sure just single nucleotide to single nucleotide
  print(nrow(hb_all_snp))
  set.seed(nsnp)
  if (nrow(hb_all_snp)>nsnp){
  hb<-hb_all_snp %>% sample_n(nsnp)
  } else {
  hb<-hb_all_snp %>% sample_n(nsnp2)
  }
  # Save list of subsetted snps for this block
  data.table::fwrite(hb, paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb",bnum,"_snps.txt"), row.names = F, quote = F, sep = "\t")
}

#### Extract the genotypes as a phased vcf file at those n snps for just white uk biobank individuals, coded as 0/1
cd /oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data
ml biology
ml plink/2.0a2
plink2 --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --export vcf --extract hb1_snps.txt --output-missing-genotype 2 --out hb1_genotypes01

plink2 --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --export vcf --extract hb2_snps.txt --output-missing-genotype 2 --out hb2_genotypes01

plink2 --pgen /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pgen \
  --pvar /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.pvar \
  --psam /oak/stanford/groups/mrivas/ukbb24983/cal/phased/ukb24983_hap_chr6_v2.psam \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --export vcf --extract hb3_snps.txt --output-missing-genotype 2 --out hb3_genotypes01

#### Process the vcf files to extract the phased data per person into haplotypes
############# MUNGE VCF
ml R/4.0
R

library(dplyr)
blocks=1:3
nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"

for (bnum in blocks){
v <- data.table::fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb",bnum,"_genotypes01.vcf"))
data.table::fwrite(v, paste0(nsnp_dir,"hb",bnum,"_vcftodf.txt"), row.names = F, quote = F, sep = "\t")
}

############# SPLIT PHASING INTO SEPARATE COLUMNS
cd /oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/

### BLOCK 1
# Define files
input="hb1_vcftodf.txt"
output="hb1_vcftodf_phased.txt"

# Generate header
head -n 1 $input | tr '\t' '\n' | awk '{print (NR<10)?$0:$0"_a\t"$0"_b"}' | tr '\n' '\t' | sed 's/\t$/\n/' > $output

# Split phasing
awk 'BEGIN{ FS=OFS="\t" } {for (i=10; i<=NF; i++) {split($i, a, "|"); $i=a[1] OFS a[2]}} 1' $input  | tail -n +2 >> $output

### BLOCK 2
# Define files
input="hb2_vcftodf.txt"
output="hb2_vcftodf_phased.txt"

# Generate header
head -n 1 $input | tr '\t' '\n' | awk '{print (NR<10)?$0:$0"_a\t"$0"_b"}' | tr '\n' '\t' | sed 's/\t$/\n/' > $output

# Split phasing
awk 'BEGIN{ FS=OFS="\t" } {for (i=10; i<=NF; i++) {split($i, a, "|"); $i=a[1] OFS a[2]}} 1' $input  | tail -n +2 >> $output

### BLOCK 3
# Define files
input="hb3_vcftodf.txt"
output="hb3_vcftodf_phased.txt"

# Generate header
head -n 1 $input | tr '\t' '\n' | awk '{print (NR<10)?$0:$0"_a\t"$0"_b"}' | tr '\n' '\t' | sed 's/\t$/\n/' > $output

# Split phasing
awk 'BEGIN{ FS=OFS="\t" } {for (i=10; i<=NF; i++) {split($i, a, "|"); $i=a[1] OFS a[2]}} 1' $input  | tail -n +2 >> $output

######################## Define haplotypes
ml R/4.0
R
library(dplyr)

nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"
blocks=c(3,2,1)

for (bnum in blocks){
  # load data
  v <- data.table::fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb",bnum,"_vcftodf_phased.txt"))
  
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
  
  # clear space
  v <- 1
  haplo1a <- 1
  haplo1b <- 1
  haplo <- 1
  mmh <- 1
  
  # hold out really rare haplotypes for now
  min_cutoff <- 10
  rarehap <- filter(m_full,total<=min_cutoff)
  m_full <- filter(m_full,total>min_cutoff)
  
  #### MAKE TREE
  # Calculate distance matrix on the haplotypes
  df <- m_full %>% ungroup %>% select(-c(hap,total,unippl,twodose,onedose,dendroname)) # Make matrix just of 0s and 1s by selecting just var IDs
  print(dim(df)) # check to make sure selected correct number columns
  uniqnames <- make.names(m_full$dendroname,unique=TRUE)
  rownames(df) <- gsub("Hap\\.\\.(\\d+)", "Hap: \\1", uniqnames)
  dm <- as.data.frame(as.matrix(dist(df, method = "binary"))) # Compute the dissimilarity matrix
  
  # Initialize
  dfs_tosplit <- list(dm)
  i=1
  nclust=0
  ngroups=2
  final_clusters <- list()
  mx <- max(m_full$total) # max individuals in a single haplotype
  # set max number of doses in a cluster before splitting more
  if (mx > 80000){ # if there is one haplotype that is super large, set that as the max
    max_threshold=mx
  } else{
    max_threshold=80000 # otherwise set 80000 as the max
  }
  min_threshold=1000 # min number of doses in a cluster (if will be below this then dont split!)
  
  print(paste0("I've been here 345: block",bnum))
  
  while (i <= length(dfs_tosplit)){
    # Go through one of the dataframes need to split
    full_df <- dfs_tosplit[[i]]
    
    if (nrow(full_df)<=1){ # if there is just one haplotype left, no further splitting
      nclust=nclust+1 # count how many have so far of these "final clusters"
      m_filt <- m_full %>% filter(dendroname %in% colnames(full_df)) # filter original dataframe to the hap in this cluster
      final_clusters[[nclust]] <- m_filt %>% mutate(final_cluster=nclust) # save that cluster to the final cluster list
    } else {
      hc <- hclust(as.dist(full_df), method = "complete") # Perform hierarchical clustering
      cluster_labels <- cutree(hc, k = ngroups)
      split <- data.frame(newcluster=cluster_labels,dendroname=names(cluster_labels))
      firstsplt <- filter(split,newcluster==1)
      secondsplt <- filter(split,newcluster==2)
      split1 <- filter(m_full,dendroname %in% firstsplt$dendroname) # filter og dataframe to match these ids
      split2 <- filter(m_full,dendroname %in% secondsplt$dendroname) # filter og dataframe to match these ids
      
      if (sum(split1$total) < min_threshold | sum(split2$total) < min_threshold){ # if either cluster is below min_threshold then don't split
        nclust=nclust+1 # count how many have so far of these "final clusters"
        m_filt <- m_full %>% filter(dendroname %in% colnames(full_df)) # filter original dataframe to the hap in this cluster
        final_clusters[[nclust]] <- m_filt %>% mutate(final_cluster=nclust) # save that cluster to the final cluster list
        
      } else {
        # Check both new dfs
        split_data <- list(split1,split2)
        for (split in split_data){
          if (sum(split$total)<max_threshold){ # if that df has total doses < max_threshold then it's done splitting
            nclust=nclust+1 # count how many have so far of these "final clusters"
            final_clusters[[nclust]] <- split %>% mutate(final_cluster=nclust) # save that cluster to the final cluster list
          } else { # otherwise keep splitting
            dsum <- filter(full_df %>% select(split$dendroname),row.names(full_df) %in% split$dendroname) # filter dist matrix to relevant haps
            dfs_tosplit <- append(dfs_tosplit,list(dsum)) # otherwise add it to the list of dfs to split
          }
          #}
        }}
    }
    i=i+1 # increase total number of times looped through
  }
  
  all <- data.table::rbindlist(final_clusters) # combine all
  all %>% group_by(final_cluster) %>% summarize(sum(total))
  
  #### Add back in the rare haps
  # Identify the most common haplotype in each cluster
  nkeep=1
  subset_df <- ungroup(all %>% group_by(final_cluster) %>% arrange(-total) %>% slice(1:nkeep))
  sh <- subset_df %>% select(-c(hap,total,unippl,twodose,onedose,dendroname,final_cluster))
  
  # Get list of the haplotypes that were dropped
  rh <- rarehap %>% select(-c(hap,total,unippl,twodose,onedose,dendroname))
  clustlist <- list()
  
  # Calculate nearest cluster
  for (i in 1:nrow(rh)){
    combh <- bind_rows(rh[i,],sh)
    cd <- as.matrix(dist(combh,method=c("binary")))
    cd[cd==0] <- 100000 # want to ignore all rows with themselves
    rnum <- as.matrix(apply(cd,1,which.min))[1]
    clustlist <- append(clustlist,rnum)
  }
  rarehap$final_cluster <- unlist(clustlist) # get cluster assignments for rare haplotypes
  all <- bind_rows(rarehap,all) # add back in the rare haplotypes
  all %>% group_by(final_cluster) %>% summarize(sum(total))
  
  # Clear up space
  dfs_tosplit <- 1
  dsum <- 1
  final_clusters <- 1
  dm <- 1
  dis <- 1
  df <- 1
  m_full <- 1
  
  # saved df with munged genotypes
  data.table::fwrite(all, paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01.txt"), row.names = F, quote = F, sep = "\t")
  
  ### MUNGE FORMAT FOR REGRESSION
  #all <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01.txt")) # read in the above if not running script from start
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
  data.table::fwrite(cassign, paste0(nsnp_dir,"hb",bnum,"_clusterdoses.txt"), row.names = F, quote = F, sep = "\t")
  cassign <- cassign %>% select(-UKB_ID)
  
  # Calculate doses
  cdoses <- data.frame(
    Cluster = 1:ncol(cassign),
    cluster_onedose = colSums(apply(cassign, 2, function(x) x == 1)),
    cluster_twodose = colSums(apply(cassign, 2, function(x) x == 2)),
    cluster_unippl = colSums(apply(cassign, 2, function(x) x %in% c(1, 2)))
  ) %>% mutate(cluster_total=cluster_twodose*2+cluster_onedose)
  
  print(paste0("I've been here 678: block",bnum))
  data.table::fwrite(cdoses, paste0(nsnp_dir,"hb",bnum,"_clusterstats.txt"), row.names = F, quote = F, sep = "\t")
  data.table::fwrite(l, paste0(nsnp_dir,"hb",bnum,"_haplostats.txt"), row.names = F, quote = F, sep = "\t")
  
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
  
  data.table::fwrite(hap_an, paste0(nsnp_dir,"hb",bnum,"_temphapan.txt"), row.names = F, quote = F, sep = "\t")
  
  # combine
  cov <- bind_cols(cov_cov,cov_traits) # combine covariate and disease info
  dat <- left_join(hap_an,cov, by=c("UKB_ID"="ID"))
  
  data.table::fwrite(dat, paste0(nsnp_dir,"hb",bnum,"_beforeregandleaveoutalldiseases.txt"), row.names = F, quote = F, sep = "\t")
  
  #### Prepare to run regressions
  # Leave out the most frequent cluster
  leaveout = paste0("cluster",(cdoses %>% filter(cluster_total==max(cdoses$cluster_total)))$Cluster)
  print(leaveout)
  keptclusters <- colnames(dat %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b,-leaveout))
  nclust <- length(keptclusters)
  
  data.table::fwrite(dat %>% select(-leaveout), paste0(nsnp_dir,"hb",bnum,"_preregressiondf_wleaveout.txt"), row.names = F, quote = F, sep = "\t") # save df before running regression
  
  # initialize data
  results_dt=data.table::data.table(trait=character(),beta=numeric(),se=numeric(),p=numeric())
  
  for (trait in colnames(cov_traits)){ # trait=mycols_traits[12]
    cov_traits %>% group_by(!!sym(trait)) %>% count()
    form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+sex+age+PC1+PC2+PC3+
      PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=dat)
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
  
  data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_afterregalldiseases.txt"), row.names = F, quote = F, sep = "\t")
  
  print(paste0("Done regression: block",bnum))

}

## Explore the data
(filter(results_dt,p<1e-6)) %>% group_by(trait) %>% count() %>% filter(n>1) 
sigtraits <- ((filter(results_dt,p<1e-6)) %>% group_by(trait) %>% count() %>% filter(n>1) )$trait
ii <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe.info.tsv")
filter(ii,`#GBE_ID` %in% sigtraits) %>% select(`#GBE_ID`,`GBE_NAME`,N)

###############################
# Plot key traits
library(dplyr)
library("dendextend")
library("RColorBrewer")
library("corrplot")
library(ComplexHeatmap)
library("ggplot2")

nsnp_dir="/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/"
z_cut <- 4
min_threshold <- 20000
blocks=1:3

for (bnum in blocks){
  ii <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.phe.info.tsv")
  r1 <- data.table::fread(paste0(nsnp_dir,"hb1_afterregalldiseases.txt")) %>% mutate(z=beta/se)
  r1sig <- unique((r1 %>% filter(abs(z)>z_cut))$trait)
  r2 <- data.table::fread(paste0(nsnp_dir,"hb2_afterregalldiseases.txt")) %>% mutate(z=beta/se)
  r2sig <- unique((r2 %>% filter(abs(z)>z_cut))$trait)
  r3 <- data.table::fread(paste0(nsnp_dir,"hb3_afterregalldiseases.txt")) %>% mutate(z=beta/se)
  r3sig <- unique((r3 %>% filter(abs(z)>z_cut))$trait)
  allsigtraits <- unique(c(r1sig,r2sig,r3sig))
  sigtraitinfo <- filter(ii,`#GBE_ID` %in% allsigtraits) %>% select(`#GBE_ID`,GBE_NAME,N) %>% rename(trait=`#GBE_ID`)
  
  # switch to wide
  if (bnum ==1) {
    h <- r1 %>% filter(trait %in% allsigtraits)
  } else if (bnum ==2){
    h <- r2 %>% filter(trait %in% allsigtraits)
  } else if (bnum == 3){
    h <- r3 %>% filter(trait %in% allsigtraits)
  }
  
  z <- data.frame(tidyr::spread(h %>% select(trait,cluster,z),cluster,z))
  traits <- z$trait
  rownames(z) <- traits
  z <- z %>% select(-trait)
  
  # filter to a subset of clusters and traits that have total doses > min_cutoff or at least 1 trait with abs(z)>z_cut
  c <- data.table::fread(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/genotype_data/hb",bnum,"_clusterstats.txt")) # cluster info, generated by hapclusterandplot.R c <- cdoses
  mckeep <- paste0("X",(c %>% filter(cluster_total>min_threshold&cluster_total!=max(c$cluster_total)))$Cluster)
  clustkeep <- unique(c(mckeep,colnames(z[,(colSums(abs(z)>z_cut))>0])))
  z_filt <- z[(rowSums(abs(z)>z_cut))>0,] # filter to a small number of traits of interest (have >= 1 clusters w/ abs(z)>z_cut in the clusters to plot)
  z_filt <- as.data.frame(z_filt) %>% select(all_of(clustkeep))
  
  #### Plot haplotypes and get cluster order for heatmap
  nkeep=40 # select the n most common haplotypes from each cluster, only clusters to keep
  m <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01.txt")) # load in genptyoe info, generated by hapclusterandplot.R
  subset_df <- m %>% group_by(final_cluster) %>% arrange(-total) %>% slice(1:nkeep)
  subset_df <- subset_df %>% filter(final_cluster %in% gsub("X","",clustkeep))
  df <- t(as.matrix(subset_df %>% ungroup %>% select(-c(hap,total,unippl,twodose,onedose,dendroname,final_cluster))))
  l<- as.character(subset_df$final_cluster)
  C <- cluster_within_group(df,l)
  set.seed(123)
  tiff(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/figures/hb",bnum,"hapdendro.tiff"),height=10,width=10,res=300,units="in")
  Heatmap(t(df), cluster_rows = C, border=TRUE, row_split = length(unique(subset_df$final_cluster)),
          show_column_dend = FALSE,
          cluster_columns = FALSE,show_column_names=FALSE,
          show_row_names = TRUE,
          col=c("black","white"),show_heatmap_legend = FALSE,
          row_gap=unit(0.8,"mm"),
          row_title_gp = gpar(fontsize = 25),
          row_title_rot = 0)#,left_annotation = rowAnnotation(Cluster = l,show_legend=F))
  dev.off()
  clustorder <- unique(l[order.dendrogram(C)]) # top to bottom
  
  #### Plot associations if hap group N > 20,000 or p < 1e-6
  # Reorder and name the clusters and subset down just to those to plot
  z <- z_filt
  z_filt <- 1
  clusternames <- paste0("X",rev(clustorder))
  z <- as.data.frame(z) %>% select(clusternames) # reorder columns
  znocut <- z
  
  # Add cutoff values for plotting
  z[z<(-5)] <- -5
  z[z>(5)] <- 5
  
  # Cluster the traits axis of the heatmap (keep clusters in order of dendrogram not clustered here
  hc <- hclust(dist(z,method="euclidean"),method="ward.D")
  z <- z[hc$order,]
  z$traits <- rownames(z)
  z$traits <- factor(z$traits,levels=traits)
  melted_df <- reshape2::melt(z)
  
  # Remove redundant traits
  melted_df <- melted_df %>% filter(!(traits %in% c("BIN_FC3006144","BIN_FC50001468","BIN_FC40001448")))
  
  # Add in updated trait names
  melted_df <- melted_df %>% left_join(sigtraitinfo,by=c("traits"="trait"))
  melted_df$GBE_NAME <- factor(melted_df$GBE_NAME,levels=unique(melted_df$GBE_NAME))
  trait_annotations <- data.frame(traits=c("HC703", "BIN_FC10006152", "BIN_FC5006152", "HC310", "HC990", "BIN_FC20020549", "HC1345", "HC151", "HC810", "HC528", "HC530", "HC170", "HC639", "BIN22133", "HC1201", "BIN21068", "BIN_FC6006149", "HC55", "HC645", "HC1132", "HC303", "HC1212", "HC430", "HC1242", "HC422", "BIN_FC9006154", "BIN_FC2006154", "BIN_FC8006154", "HC1188", "BIN_FC1006149", "HC1102", "HC201", "BIN10030820", "HC38", "BIN_FC10002986", "HC219", "HC643","HC1021", "HC49", "HC85", "FH1044", "HC1094", "HC1322", "BIN_FC9006152", "HC1581", "HC382", "BIN_FC4006152", "HC1036", "BIN_FC3006153", "BIN_FC8006153", "HC648", "BIN_FC3006177", "BIN_FC6006177", "FH1220", "BIN2443", "HC221", "HC652", "BIN23067", "BIN2188", "BIN_FC4001747", "HC1233", "HC352", "HC125", "HC607", "HC497", "BIN2473", "HC620", "HC868", "HC1209", "HC26", "cancer1003", "cancer1060", "HC646", "HC332", "BIN_FC10001707", "HC1159", "HC1213", "HC91","HC321", "HC185", "HC1158", "BIN23069", "HC900", "cancer1041", "HC853", "BIN22127", "BIN2453", "BIN1950", "FH1065", "HC215", "HC700", "BIN_FC1006153", "BIN_FC6006153", "BIN23066", "HC337", "HC96", "HC156", "HC1024", "HC78", "HC136", "BIN_FC10002492", "HC1236", "BIN2316", "HC1211"),Plot=c("Mineral metabolism disorder","Rhinitis/eczema","Rhinitis/eczema 2","Hiatal hernia","Varicose veins","Anxiety Rx","Vulvovaginitis","MS","MS 2","HSV infection","VZV infection","Sarcoidosis","Sarcoidosis 2","Sarcoidosis 3","SLE","Celiac disease","Dentures","Hyperthyroid thyrotoxicosis","Hyperthyroid thyrotoxicosis 2","Malabsorption","Celiac disease 2","Rheumatoid arthritis 2","Rheumatoid arthritis","Ankylosing spondylitis 2","Ankylosing spondylitis", "Paracetamol use","Ibuprofen use","Ibuprofen use 2","Vitiligo","Mouth ulcers","UC 2","UC","Elevated rheumatoid factor","Psoriasis","Insulin","Hypothyroidism","Hypothyroidism 2","Allergic rhinitis","Allergic rhinitis 2","Enlarged prostate","Prostate cancer","Inguinal hernia","Enlarged prostate 2","Asthma","Asthma 2","Asthma 3","Asthma 4","Asthma 5","Diabetes Rx (F)","Insulin 2","Diabetes 5","Diabetes Rx (M)","Insulin 3","Diabetes 1","Diabetes 2","Diabetes 3","Diabetes 4","MCV seropositivity","Chronic illness","Brown hair","SLE 2","SLE 3","Anemia","Anemia 2","STI","Other serious disability","Anemia 3","Iridocyclitis","Reactive arthropathies","Iritis","Skin cancer","Non-melanoma skin cancer","Thyroiditis","PTSD","Right handedness","Psoriasis 2","Psoriatic arthropathies","Psoriatic arthropathy 2","Sjogren's","Cellulitis","Dermatitis","HPV18 seropositivity","Otitis externa","Cervical cancer","Stye","Asthma 6","Cancer","Sensitivity","Hypertension 2","Hypertension","Lipidemia","Cholesterol Rx (F)","Cholesterol Rx","JCV seropositivity","Type 1 Diabetes","Polymyalgia rheumatica","Nasal polyps","Nasal polyps 2","Arthritis","Connective tissue disorder","Other Rx","Connective tissue disorder 2","Wheezing","Seropositive RA"))
  melted_df <- left_join(melted_df,trait_annotations,by=c("traits"))
  
  # Remove redundant traits and rename remaining accordingly
  if (bnum == 1) {
    melted_df <- melted_df %>% filter(!(traits %in% c("HC170","HC639","HC1132","BIN_FC8006154","BIN_FC2006154")))
    melted_df$Plot[melted_df$traits == "BIN22133"] <- "Sarcoidosis"
    melted_df$Plot[melted_df$traits == "HC303"] <- "Malabsorption"
    melted_df <- melted_df[!grepl(" (2|3)$", melted_df$Plot), ]
  } else if (bnum == 2){
    melted_df <- melted_df %>% filter(!(traits %in% c("BIN_FC3006153","BIN_FC8006153","BIN_FC10002986","HC648","BIN_FC3006177",
                                                      "FH1220","HC221","HC652","HC1132","HC170","HC639","BIN_FC10001707","BIN_FC4001747","HC1233")))
    melted_df$Plot[melted_df$traits == "BIN22133"] <- "Sarcoidosis"
    melted_df$Plot[melted_df$traits == "HC303"] <- "Malabsorption"
    melted_df$Plot[melted_df$traits == "BIN2443"] <- "Diabetes"
    melted_df$Plot[melted_df$traits == "HC352"] <- "SLE"
    melted_df <- melted_df[!grepl(" (2|3|4|5)$", melted_df$Plot), ]
  } else if (bnum == 3){
    melted_df <- melted_df %>% filter(!(traits %in% c("HC1021","HC1132","HC170","HC639","BIN_FC9006152","BIN1950","BIN_FC1006153","HC700","BIN_FC3006153", "BIN_FC8006153","BIN_FC10002986", "HC648",
                "BIN_FC3006177","FH1220", "HC221","HC652","HC337")))
    melted_df$Plot[melted_df$traits == "BIN22133"] <- "Sarcoidosis"
    melted_df$Plot[melted_df$traits == "HC303"] <- "Malabsorption"
    melted_df$Plot[melted_df$traits == "BIN2443"] <- "Diabetes"
    melted_df$Plot[melted_df$traits == "HC49"] <- "Allergic rhinitis"
    melted_df$Plot[melted_df$traits == "BIN22127"] <- "Asthma"
    melted_df$Plot[melted_df$traits == "BIN_FC6006153"] <- "Cholesterol lowering Rx"
    melted_df <- melted_df[!grepl(" (2|3|4|5)$", melted_df$Plot), ]
  }
  
  # Set factors
  melted_df$Plot <- factor(melted_df$Plot,levels=unique(melted_df$Plot))
  
  # Reformat hap groups
  melted_df$variable <- gsub("X","",melted_df$variable)
  melted_df$variable <- factor(melted_df$variable,levels=unique(melted_df$variable))
  
  # Plot the heatmap using ggplot2
  ggplot(melted_df, aes(x = Plot, y = variable, fill = value)) + # x=traits
    geom_tile()+
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits = c(-5, 5))+
    theme_classic()+
    theme(text = element_text(size=25),
          axis.text.x = element_text(angle = -70, hjust = 0))+
    #labs(x="Traits",y="Clusters",fill="Z-score")
    labs(x="",y="",fill="Z-score")+
    scale_x_discrete(limits=unique(melted_df$Plot)) # z$traits
  ggsave(paste0("/oak/stanford/groups/pritch/users/courtrun/hla_ukbb/figures/hb",bnum,"hapheatmap.tiff"),height=10,width=20)
  
  # Save plotted data
  data.table::fwrite(melted_df %>% rename(hapgroup=variable,Z=value,LONG_NAME=GBE_NAME) %>% mutate(block=bnum),paste0(nsnp_dir,"hb",bnum,"_heatmapplotted.txt"), row.names = F, quote = F, sep = ",")
  
  # Save uncut data
  znocut$traits <- rownames(znocut)
  melted_df_nocut <- reshape2::melt(znocut)
  melted_df_nocut <- melted_df_nocut %>% filter(traits %in% melted_df$traits) # same traits as in plot
  melted_df_nocut <- left_join(melted_df_nocut,trait_annotations,by=c("traits")) # add annotation info
  melted_df_nocut <- melted_df_nocut %>% left_join(sigtraitinfo,by=c("traits"="trait")) # add annotation info
  melted_df_nocut$variable <- gsub("X","",melted_df_nocut$variable) # reformat hap groups
  
  data.table::fwrite(melted_df_nocut %>% rename(hapgroup=variable,Z=value,LONG_NAME=GBE_NAME) %>%
                       mutate(block=bnum) %>%
                       select(traits,hapgroup,Z,LONG_NAME,N,Plot,block),paste0(nsnp_dir,"hb",bnum,"_heatmapplotted_nocut.txt"), row.names = F, quote = F, sep = ",")
  
}


