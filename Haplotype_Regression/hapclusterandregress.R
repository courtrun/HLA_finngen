## Written by Courtney Smith 7-1-2023
## Goal of script: Defining clusters and running haplotype regression pipeline - Blocks 1-3 - 30,100,150 SNPs

library(dplyr)

an <- data.table::fread("/home/ivm/from_satu/general_files/rsids.txt.gz")

blocks=c("1","2","3")

for (bnum in blocks) {
  snps <- data.table::fread(paste0("/home/ivm/from_satu/hbg_full_sets/hbg",bnum,"_full.tsv"),header=F)
  filt <- filter(an,rsids %in% snps$V1 & rsids!="")
  filtan <- filt %>% mutate(V1=paste0("chr",`#chrom`,"_",pos,"_",ref,"_",alt)) %>% select(V1)
  data.table::fwrite(filtan, paste0("/home/ivm/from_satu/hbg_full_sets/hbg",bnum,"_maf.txt"), row.names = F, col.names=F, quote = F, sep = "\t")
}

##
library("dplyr")
library("dendextend")
library("RColorBrewer")
library(ggplot2)
library(ComplexHeatmap)
library(cowplot)
library(stringdist)
library(data.table) # note: this has a "set" function that conflicts with dendextend

# Set parameters
blocks = c(1,2,3,"2b","3b")
nsnps= c(100,1000) # c(30,100,120,130,140,150)
############ Select subset of n random snps from the full set of snps in the block

for (nsnp in nsnps){
  for (bnum in blocks){
    hla_dir="/home/ivm/haplo_analysis/"
    nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
    setwd(nsnp_dir)

    # Load in snps in the block, randomly subset
    hb_all<-fread(paste0("/home/ivm/from_satu/hbg_full_sets/hbg",bnum,"_maf.txt"), header = F) # all snps in the block w/ MAF > 1%
    hb_all_snp <- hb_all[grep("^chr6_\\d+_[A-Z]_[A-Z]$", hb_all$V1), , drop = FALSE] # make sure just single nucleotide to single nucleotide
    set.seed(nsnp)
    hb<-hb_all_snp %>% sample_n(nsnp)
    #hb<-hb_all %>% sample_n(nsnp)

    # Save list of subsetted snps for this block
    data.table::fwrite(hb, paste0(nsnp_dir,"hb",bnum,"_snps.txt"), row.names = F, quote = F, sep = "\t")
  }
}

########### Code for extracting the haplotype individual data for these snps on the terminal in dif formats
####### Create the genotype file for the snps in this block with the alleles (A,T,C,Gs)
##extract the snps with PLINK
##cd /home/ivm/haplo_analysis/snps100
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb3_snps.txt --out hb3_genotypes
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb2_snps.txt --out hb2_genotypes
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb1_snps.txt --out hb1_genotypes

##cd /home/ivm/haplo_analysis/snps30
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb3_snps.txt --out hb3_genotypes
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb2_snps.txt --out hb2_genotypes
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb1_snps.txt --out hb1_genotypes

##cd /home/ivm/haplo_analysis/snps120
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb3_snps.txt --out hb3_genotypes
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb2_snps.txt --out hb2_genotypes
##plink --bfile /finngen/library-red/finngen_R10/genotype_plink_1.0/data/finngen_R10 --keep-allele-order --keep /home/ivm/from_satu/general_files/individuals2.txt --recode --extract hb1_snps.txt --out hb1_genotypes

###### Create the genotype file for the snps in this block as 0s and 1s
##cd /home/ivm/haplo_analysis/snps30 and repeat for snps100 and for snps120
##plink --recode 01 transpose --file hb1_genotypes --output-missing-genotype 2 --out hb1_genotypes01
##plink --recode 01 transpose --file hb2_genotypes --output-missing-genotype 2 --out hb2_genotypes01
##plink --recode 01 transpose --file hb3_genotypes --output-missing-genotype 2 --out hb3_genotypes01
##plink --recode 01 --file hb1_genotypes --output-missing-genotype 2 --out hb1_genotypes01
##plink --recode 01 --file hb2_genotypes --output-missing-genotype 2 --out hb2_genotypes01
##plink --recode 01 --file hb3_genotypes --output-missing-genotype 2 --out hb3_genotypes01

########### Define the clusters and run regression

library(dplyr)
#library(stringdist)
library(data.table)
library("RColorBrewer")

# Set parameters
blocks = c(1,2,3,"2b","3b")
nsnps= c(100,1000) # c(100,500,1000) # c(30,100,120,130)

for (nsnp in nsnps){
  for (bnum in blocks){
    hla_dir="/home/ivm/haplo_analysis/"
    nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
    setwd(nsnp_dir)

    print(paste0("Just started: block",bnum," and nsnp",nsnp))

    ## Load data and munge to phased
    indiv_data <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_genotypes01.tped"))
    temp <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_genotypes01.ped"))
    finngenid <- temp$V2

    selected_cols <- c(2,seq(5, ncol(indiv_data), by = 2))
    haplo1a <- indiv_data[, ..selected_cols]
    rownames(haplo1a) <- haplo1a$V2
    haplo1a <- haplo1a %>% select(-V2)
    snps <- rownames(haplo1a)
    haplo1at <- as.data.frame(t(haplo1a))
    colnames(haplo1at) <- snps

    selected_cols <- c(2,seq(6, ncol(indiv_data), by = 2))
    haplo1b <- indiv_data[, ..selected_cols]
    rownames(haplo1b) <- haplo1b$V2
    haplo1b <- haplo1b %>% select(-V2)
    snps <- rownames(haplo1b)
    haplo1bt <- as.data.frame(t(haplo1b))
    colnames(haplo1bt) <- snps
    haplo <- bind_rows(haplo1at,haplo1bt)
    m <- distinct(haplo) # df with n snps as columns and 2x individuals as columns, keeping only distinct haplotypes

    # Add rownames with haplo info
    mh <- tidyr::unite(m,hap,everything(),sep="") # create a new column "hap" that combines all other columns
    rownames(m) <- mh$hap

    # Make a dataframe that is m with extra columns like hap numbers
    # Create dataframes with info in alt ways
    h <- tidyr::unite(haplo,hap,everything(),sep="") # create a new column "hap" that combines all other columns
    first <- data.frame(ID=paste0("ID ",c(1:(nrow(h)/2))),haplo1a=(h[c(1:(nrow(h)/2)),]),haplo1b=(h[c((nrow(h)/2+1):(nrow(h))),]))
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

    print(paste0("I've been here 123: block",bnum," and nsnp",nsnp))

    # Open space
    indiv_data <- 1
    temp <- 1
    haplo1a <- 1
    haplo1b <- 1
    haplo <- 1
    mmh <- 1

    # Only keep the hap that have > min_cutoff <- uncomment this section if want to filt to reasonable dist matrix for lots of snps
    if (nsnp>100){
      min_cutoff <- 10
      rarehap <- filter(m_full,total<=min_cutoff)
      m_full <- filter(m_full,total>min_cutoff)
    }

    # Calculate distance matrix on full data
    df <- m_full %>% ungroup %>% select(matches("_")) # Make matrix just of 0s and 1s
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
    #min_threshold=0 # min number of doses in a cluster <- for now not including minimum

    print(paste0("I've been here 345: block",bnum," and nsnp",nsnp))

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
        #if (sum(split1$total) < min_threshold | sum(split2$total) < min_threshold){ # if either cluster is below min_threshold then don't split
        # nclust=nclust+1 # count how many have so far of these "final clusters"
        #m_filt <- m_full %>% filter(dendroname %in% colnames(full_df)) # filter original dataframe to the hap in this cluster
        #final_clusters[[nclust]] <- m_filt %>% mutate(final_cluster=nclust) # save that cluster to the final cluster list
        #} else {
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
        }
      }
      i=i+1 # increase total number of times looped through
    }

    all <- rbindlist(final_clusters) # combine all
    #tail(all %>% select(hap,total,final_cluster),10)
    #dim(all %>% group_by(final_cluster) %>% summarize(sum(total))) # how many clusters
    #head(all %>% group_by(final_cluster) %>% summarize(sum=sum(total)) %>% arrange(-sum))

    ### Add back in the rare haps
    # Identify the most common haplotype in each cluster
    nkeep=1
    subset_df <- ungroup(all %>% group_by(final_cluster) %>% arrange(-total) %>% slice(1:nkeep))
    sh <- subset_df %>% select(matches("_"),-final_cluster)

    # Get list of the haplotypes that were dropped
    rh <- rarehap %>% select(matches("_"))
    clustlist <- list()

    # Calculate nearest cluster
    for (i in 1:nrow(rh)){
      combh <- bind_rows(rh[i,],sh)
      cd <- as.matrix(dist(combh,method=c("binary")))
      cd[cd==0] <- 100000 # want to ignore all rows with themselves
      rnum <- as.matrix(apply(cd,1,which.min))[1]
      clustlist <- append(clustlist,rnum)
    }

    rarehap$final_cluster <- unlist(clustlist)

    all <- bind_rows(rarehap,all)

    # Clear up space
    dfs_tosplit <- 1
    dsum <- 1
    final_clusters <- 1
    dm <- 1
    dis <- 1
    df <- 1
    m_full <- 1

    print(paste0("I've been here 456: block",bnum," and nsnp",nsnp))

    # Add color scheme for all clusters with more than x haplotypes
    allplot <- all %>% group_by(final_cluster) %>% summarize(sum=sum(total)) %>% arrange(-sum) %>% filter(sum>5000)
    nplot <- nrow(allplot) # how many clusters have above 10,000 doses

    # Colors
    #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector = c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00","black", "gold1", "skyblue2", "#FB9A99", "palegreen2",  "#CAB2D6", "#FDBF6F", "gray70", "khaki2",  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
    colors <- data.frame(final_cluster=unique(allplot$final_cluster),
                         color= sample(col_vector,size=nplot))

    # Reorder to dendrogram and get colors
    all <- left_join(all,colors,by=c("final_cluster"))
    all$dendroname <- paste0("Cluster ",all$final_cluster,": ",all$total)

    # saved df with munged genotypes
    data.table::fwrite(all, paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01.txt"), row.names = F, quote = F, sep = "\t")

    ### MUNGE FORMAT FOR REGRESSION
    l <- all %>% rename(Cluster=final_cluster,Haplotype=hap) # long form of haplotypes, num individudals with that haplotype, cluster that haplotype belongs to
    temp<-left_join(first,l,by=c("haplo1a"="Haplotype")) %>% rename(cluster1a=Cluster)
    hap_an<-left_join(temp, l, by=c("haplo1b"="Haplotype")) %>% rename(cluster1b=Cluster)
    hap_an$FINNGENID <- finngenid
    hap_an <- filter(hap_an,!is.na(cluster1b)&!is.na(cluster1a))
    print(dim(hap_an))
    nclust=max(hap_an$cluster1a,hap_an$cluster1b)
    for (n in 1:nclust){
      varname = paste0("cluster",n)
      hap_an[, varname] <- ifelse(hap_an$cluster1a==n&hap_an$cluster1b==n,2,ifelse(hap_an$cluster1a==n|hap_an$cluster1b==n,1,0))
    }

    print(paste0("I've been here 567: block",bnum," and nsnp",nsnp))

    # Make df of just ID and cluster belonging
    cassign <- hap_an %>% select(FINNGENID,matches("cluster")) %>% select(-cluster1a,-cluster1b)
    data.table::fwrite(cassign, paste0(nsnp_dir,"hb",bnum,"_clusterdoses.txt"), row.names = F, quote = F, sep = "\t")

    cassign <- cassign %>% select(-FINNGENID)

    # Calculate doses
    cdoses <- data.frame(
      Cluster = 1:ncol(cassign),
      cluster_onedose = colSums(apply(cassign, 2, function(x) x == 1)),
      cluster_twodose = colSums(apply(cassign, 2, function(x) x == 2)),
      cluster_unippl = colSums(apply(cassign, 2, function(x) x %in% c(1, 2)))
    ) %>% mutate(cluster_total=cluster_twodose*2+cluster_onedose)

    print(paste0("I've been here 678: block",bnum," and nsnp",nsnp))

    # Add in covariate info
    #covariates and traits
    cov<-fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_PHENO_V1.txt.gz")
    traits<-fread("/home/ivm/from_satu/general_files/R10_condit_results.txt")
    traits<-traits[!duplicated(traits$trait),]
    traits<-traits[,c(trait)]
    rest<-c("FINNGENID","AGE_AT_DEATH_OR_END_OF_FOLLOWUP","SEX_IMPUTED","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
    vars<-append(rest, traits)
    dat<-cov[,vars, with=F]
    dat[is.na(dat)]<-0
    dat<-merge(hap_an,dat, by.y = "FINNGENID", by.x = "FINNGENID", all.y = T)

    ### Save files have so far before running regression
    data.table::fwrite(cdoses, paste0(nsnp_dir,"hb",bnum,"_clusterstats.txt"), row.names = F, quote = F, sep = "\t")
    data.table::fwrite(l, paste0(nsnp_dir,"hb",bnum,"_haplostats.txt"), row.names = F, quote = F, sep = "\t")
    #data.table::fwrite(dat, paste0(nsnp_dir,"hb",bnum,"_fullregdf.txt"), row.names = F, quote = F, sep = "\t")

    # Leave out the most frequent cluster
    leaveout = paste0("cluster",(cdoses %>% filter(cluster_total==max(cdoses$cluster_total)))$Cluster)
    print(leaveout)
    keptclusters <- colnames(dat %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b,-leaveout))
    nclust <- length(keptclusters)

    ## Cluster info
    ninfo <- length(colnames(hap_an))+length(rest) - 1

    ## binary traits
    dat_bin <- dat %>% select(-c("BMI_IRN","WEIGHT_IRN","HEIGHT_IRN"))

    ## continious
    dat_con <- dat %>% select(c(1:ninfo),"BMI_IRN","WEIGHT_IRN","HEIGHT_IRN")

    print(paste0("I've been here 789: block",bnum," and nsnp",nsnp))

    data.table::fwrite(dat_bin %>% select(-leaveout), paste0(nsnp_dir,"hb",bnum,"_preregressiondf.txt"), row.names = F, quote = F, sep = "\t") # save df before running regression

    ##regression analysis, remember to exclude the most common cluster!
    ## for continuous traits the regression can be run similarly, but family="binomial" needs to be taken out

    results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())
    bin_traits <- colnames(dat_bin %>% select(-c(1:ninfo)))

    print(paste0("About to start regression: block",bnum," and nsnp",nsnp))

    # Run regression
    for (trait in bin_traits) {
      #log<-glm(T1D~SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+cluster1+cluster3+cluster4+cluster5+cluster6+cluster7+cluster8+cluster9+cluster10, family="binomial", data=dat)
      ##=summary(glm(dat[[j]]~0+cluster1+cluster2+cluster3+cluster4+cluster5+cluster6+cluster7+cluster8+cluster9+cluster10+cluster11+cluster12+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=dat))$coefficients[1:12,1]
      form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
      regres <- glm(form, family="binomial", data=dat_bin)
      beta<-summary(regres)$coefficients[2:(nclust+1),1]
      se<-summary(regres)$coefficients[2:(nclust+1),2]
      p<-summary(regres)$coefficients[2:(nclust+1),4]
      tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
      results_dt=rbindlist(list(results_dt,tmp_dt))
    }

    # save regression results
    cluster<-gsub("cluster","",keptclusters)
    results_dt$cluster<-rep(cluster, length(bin_traits))
    data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_regresults.txt"), row.names = F, quote = F, sep = "\t")

    print(paste0("Done regression: block",bnum," and nsnp",nsnp))

    ##continious traits regression
    #results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())
    #cont_traits <- colnames(dat_con %>% select(-c(1:ninfo)))

    #for (trait in cont_traits) {
    ##log<-glm(T1D~SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+cluster1+cluster3+cluster4+cluster5+cluster6+cluster7+cluster8+cluster9+cluster10, family="binomial", data=dat)
    ##=summary(glm(dat[[j]]~0+cluster1+cluster2+cluster3+cluster4+cluster5+cluster6+cluster7+cluster8+cluster9+cluster10+cluster11+cluster12+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, family="binomial", data=dat))$coefficients[1:12,1]
    # form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    #regres <- glm(form, data=dat_con)
    #beta<-summary(regres)$coefficients[2:(nclust+1),1]
    #se<-summary(regres)$coefficients[2:(nclust+1),2]
    #p<-summary(regres)$coefficients[2:(nclust+1),4]
    #tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
    #results_dt=rbindlist(list(results_dt,tmp_dt))
    #}

    # save continous trait regression results
    #cluster<-gsub("cluster","",keptclusters)
    #results_dt$cluster<-rep(cluster, length(cont_traits))
    #data.table::fwrite(results_dt, paste0(nsnp_dir,"hb",bnum,"_regresultscont.txt"), row.names = F, quote = F, sep = "\t")
  }
}
