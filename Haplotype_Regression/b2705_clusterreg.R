# Courtney Smith
# Goal of script: Repeat associations analysis on same original haplotype groups but with only individuals negative for B2705

library(dplyr)
library(data.table)
library(car)
library(ggplot2)

nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

## Adjust cluster regression for an allele
blocks=c(1:3)

for (bnum in blocks){
  print(bnum)
  dat_bin <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_preregressiondf.txt"))
  ninfo <- ncol(dat_bin)-269 # only keep columns that have phenotype info, assumes 269 traits
  bin_traits <- colnames(dat_bin %>% select(-c(1:ninfo)))
  keptclusters <- colnames(dat_bin %>% select(matches("cluster")) %>% select(-cluster1a,-cluster1b))
  nclust=length(keptclusters)

  a <- data.table::fread("/finngen/red/courtrun/R10_amino_allele/R10_alleles98_MAF.txt")
  colnames(a) <- gsub("\\*|:","",colnames(a)) # reformat column names to remove special characters
  dat_bin <- left_join(dat_bin,a,by=c("FINNGENID"="FID"))

  results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric())
  dat_bin <- filter(dat_bin,B2705==0) # filter to B2705 negative
  ## ALL TRAITS
  # Run regression
  for (trait in bin_traits) {
    print(trait)
    form <- as.formula(paste(paste(trait,"~ "),paste(keptclusters,collapse="+"),paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=dat_bin)
    beta<-summary(regres)$coefficients[2:(nclust+1),1]
    se<-summary(regres)$coefficients[2:(nclust+1),2]
    p<-summary(regres)$coefficients[2:(nclust+1),4]
    tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p)
    results_dt=rbindlist(list(results_dt,tmp_dt))
    print(trait)
  }
  print("hi")
  # save regression results
  cluster<-gsub("cluster","",keptclusters)
  results_dt$cluster<-rep(cluster, length(bin_traits))
  results_dt$cluster <- as.integer(results_dt$cluster)
  data.table::fwrite(results_dt, paste0(nsnp_dir,"b2705neg_clusterreg_bnum",bnum,".txt"), row.names = F, quote = F, sep = "\t")
}
######################
# Allele regressions
nsnp=1000
hla_dir="/home/ivm/haplo_analysis/"
nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
setwd(nsnp_dir)

a <- data.table::fread("/finngen/red/courtrun/R10_amino_allele/R10_alleles98_MAF.txt")
colnames(a) <- gsub("\\*|:","",colnames(a)) # reformat column names to remove special characters

# Get traits
cov<-fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_PHENO_V1.txt.gz")
traits <- fread("/finngen/green/courtrun/finngen_R10_endpoint_core_noncore_1.0.txt") %>%
  filter(CORE_ENDPOINTS=="yes") %>% select(NAME) %>% rename(trait=NAME)
traits<-traits[!duplicated(traits$trait),]
traits<-traits[,c(trait)]
rest<-c("FINNGENID","AGE_AT_DEATH_OR_END_OF_FOLLOWUP","SEX_IMPUTED","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
vars<-append(rest, traits)
dat<-cov[,vars, with=F]
dat[is.na(dat)]<-0
bin_traits <- traits

# Df with alleles and covariates
dat<-merge(a,dat, by.y = "FINNGENID", by.x = "FID")

results_dt=data.table(nimi=character(),beta=numeric(),se=numeric(),p=numeric(),allele=character())

keptclusters <- colnames(a %>% select(-FID)) # list of alleles to run regression for
nclust <- 1

#filter to just b27 negative individuals
dat <- dat %>% filter(B2705==0)

#filter to traits sig associated with b2705
#r <- data.table::fread(paste0(nsnp_dir,"alleles_indiv_regresults_alltraits.txt"))
#b27traits <- filter(r,allele=="B2705"&p<2e-7)$nimi

for (trait in bin_traits) {
  #trait="M13_PSORIARTH"
  for (allele in keptclusters){
    form <- as.formula(paste(paste(trait,"~ "),allele,paste0("+SEX_IMPUTED+AGE_AT_DEATH_OR_END_OF_FOLLOWUP+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))
    regres <- glm(form, family="binomial", data=dat)
    beta<-summary(regres)$coefficients[2:(nclust+1),1]
    se<-summary(regres)$coefficients[2:(nclust+1),2]
    p<-summary(regres)$coefficients[2:(nclust+1),4]
    tmp_dt=data.table("nimi"=trait,"beta"=beta,"se"=se,"p"=p,allele=allele)
    results_dt=rbindlist(list(results_dt,tmp_dt))
    print(allele)
  }
  print(trait)
}

data.table::fwrite(results_dt, paste0(nsnp_dir,"alleles_indiv_regresults_alltraits_b2705neg.txt"), row.names = F, quote = F, sep = "\t")

#######################
# Analyze results

# Which traits have sig associations with alleles in b2705 negative individuals
#b <- data.table::fread(paste0(nsnp_dir,"alleles_indiv_regresults_alltraits_b2705neg.txt")) # data.table::fread(paste0(nsnp_dir,"b2705neg_clusterreg_bnum",bnum,".txt"))

# Haplotypes disease associations in b2705 negative individuals - block 1
bnum = 1
h <- data.table::fread(paste0(nsnp_dir,"b2705neg_clusterreg_bnum",bnum,".txt")) %>% rename(trait=nimi)

## Written by Courtney Smith 10-2-2024
## Goal of script: For HLA-B27 negative individuals Generate dendrogram and heatmap from (subsetted) haplotypes and heatmap of (non-subsetted) regression result

library(dplyr)
library(ggplot2)
library("dendextend")
library("RColorBrewer")
library("corrplot")
library(ComplexHeatmap)

# Set parameters
nsnps=c(1000) #c(30,100,120)
blocks=1:3
pcut <- 1e-6
min_threshold=20000
ntraits=1
zcut = 4
extended = "N"

nsnp=1000

bnum=1

for (bnum in blocks){
  #### SUBSET DOWN TO CLUSTERS AND TRAITS TO FOCUS ON
  # set directories
  hla_dir="/home/ivm/haplo_analysis/"
  nsnp_dir=nsnp_dir=paste0(hla_dir,"snps",nsnp,"/")
  setwd(nsnp_dir)

  # load data
  c <- data.table::fread(paste0("hb",bnum,"_clusterstats.txt")) # cluster info, generated by hapclusterandplot.R
  m <- data.table::fread(paste0(nsnp_dir,"hb",bnum,"_mungedgenotypes01.txt")) # load in genptyoe info, generated by hapclusterandplot.R
  h1 <- data.table::fread(paste0(nsnp_dir,"b2705neg_clusterreg_bnum",bnum,".txt")) %>% rename(trait=nimi)

  ## Get stats
  # how many unique haplotypes total
  length(unique(m$hap))

  # how many traits have at least one cluster w/ p < pcut of all clusters
  length(unique(filter(h1,p<pcut)$trait)) # (of the 269 tested)

  ## Calculate z-scores and add in dropped cluster
  h <- h1 %>% mutate(z=beta/se)
  beta <- tidyr::spread(h %>% select(trait,cluster,beta),cluster,beta)
  z <- tidyr::spread(h %>% select(trait,cluster,z),cluster,z)
  traits <- z$trait
  rownames(z) <- traits
  z <- z %>% select(-trait)
  # add in dropped cluster w/ rescale
  z <- z %>% mutate(dropped=0)
  # Function to rescale row values to have mean 0
  rescale_mean_zero <- function(row) {
    mean_zero <- row - mean(row)
    return(mean_zero)
  }
  # Apply the rescaling function to each row of the dataframe
  z_full <- t(apply(z, 1, rescale_mean_zero))
  rownames(z_full) <- traits

  # identify dropped cluster and rename last column added (dropped column) with it
  dropped <- filter(c,cluster_total==max(c$cluster_total))$Cluster
  numcol <- ncol(z_full)
  colnames(z_full) <- c(head(colnames(z_full),numcol-1),dropped)
  droppedextra <- setdiff(c$Cluster,colnames(z_full))

  # filter to a small number of traits of interest (have >= 1 clusters w/ abs(z)>zcut in the clusters to plot)
  z_filt <- z_full[(rowSums(abs(z_full)>zcut))>0,]

  # filter to a subset of clusters that have total doses > min_cutoff or at least 1 trait with abs(z)>zcut
  mckeep <- (c %>% filter(cluster_total>min_threshold))$Cluster
  mckeep <- mckeep[!mckeep %in% droppedextra]
  clustkeep <- unique(c(mckeep,colnames(z_full[,(colSums(abs(z_full)>zcut))>0])))
  z_filt <- as.data.frame(z_filt) %>% select(all_of(clustkeep))

  # stats
  dim(z_filt) # rows is traits, columns is clusters

  # save these results
  data.table::fwrite(z_filt, paste0(nsnp_dir,"hb",bnum,"_regdatatoplotb2705neg.txt"), row.names = T, quote = F, sep = "\t")
  data.table::fwrite(as.data.frame(z_full), paste0(nsnp_dir,"hb",bnum,"_regdatafullb2705neg.txt"), row.names = T, quote = F, sep = "\t")
  z_full_long <- as.data.frame(z_full)
  z_full_long$traits <- rownames(z_full)
  z_full_long <- reshape2::melt(z_full_long)
  data.table::fwrite(z_full_long %>% rename(hapgroup=variable,Z_rescaled=value), paste0(nsnp_dir,"hb",bnum,"_regdatafulllongb2705neg.txt"), row.names = T, quote = F, sep = "\t")
  ### Add color assignments
  # Step 1: Generate colors
  colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  # Step 2: Convert to RGB
  rgb_values = col2rgb(colors)
  # Step 3: Calculate brightness
  brightness = 0.299*rgb_values[1,] + 0.587*rgb_values[2,] + 0.114*rgb_values[3,]
  # Step 4: Filter colors
  threshold = 150
  col_vector = colors[brightness < threshold]
  #col_vector = c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00","black", "gold1", "skyblue2", "#FB9A99", "palegreen2",  "#CAB2D6", "#FDBF6F", "gray70", "khaki2",  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
  colors <- data.frame(cluster=colnames(z_filt),
                       color= sample(col_vector,size=ncol(z_filt)))

  #### MAKE THE SUBSETDOWN DENDROGRAM AND SNP HEATMAP
  # select the n most common haplotypes from each cluster, only clusters to keep
  nkeep=40
  subset_df <- m %>% group_by(final_cluster) %>% arrange(-total) %>% slice(1:nkeep)
  subset_df <- subset_df %>% filter(final_cluster %in% clustkeep)
  df <- t(as.matrix(subset_df %>% ungroup %>% select(matches("chr6"))))
  l<- as.character(subset_df$final_cluster)
  C <- cluster_within_group(df,l)
  set.seed(123)
  tiff(paste0(fdir,"snps",nsnp,"/figures_askout07012024/",nsnp,"_b",bnum,"minthres_",min_threshold,"_traits",ntraits,"_dendrob2705neg.tiff"),height=10,width=10,res=300,units="in")
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

  #### PLOT THE REGRESSION RESULTS HEATMAP
  # Set parameters
  fdir="/home/ivm/haplo_analysis/figures/"

  # Reorder and name the clusters and subset down just to those to plot
  z <- z_filt
  clusternames <- rev(clustorder)
  z <- as.data.frame(z) %>% select(clusternames) # reorder columns
  z_notcut <- z

  # Add cutoff values for plotting
  z[z<(-5)] <- -5
  z[z>(5)] <- 5

  # Cluster the traits axis of the heatmap (keep clusters in order of dendrogram not clustered here
  hc <- hclust(dist(z,method="euclidean"),method="ward.D")
  z <- z[hc$order,]
  z$traits <- rownames(z)
  z$traits <- factor(z$traits,levels=traits)
  melted_df <- reshape2::melt(z) # traits as id

  # Add in updated trait names
  tn <- data.table::fread("/home/ivm/from_satu/general_files/new_plot_names.txt",header=T) %>%
    select(trait,Plot)

  melted_df <-melted_df %>% left_join(tn,by=c("traits"="trait"))

  # Plot the heatmap using ggplot2
  ggplot(melted_df, aes(x = Plot, y = variable, fill = value)) + # x=traits
    #geom_tile(color="black",size=1)+
    geom_tile()+
    #scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0)+
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, limits = c(-5, 5))+
    theme_classic()+
    theme(text = element_text(size=25),
          axis.text.x = element_text(angle = -70, hjust = 0))+
    labs(x="",y="",fill="Z-score")+
    scale_x_discrete(limits=unique(melted_df$Plot)) # z$traits
  ggsave(paste0(fdir,"snps",nsnp,"/figures_askout07012024/",nsnp,"_b",bnum,"minthres_",min_threshold,"_traits",ntraits,"_cutoff5b2705neg.tiff"),height=10,width=20)
  data.table::fwrite(melted_df %>% rename(hapgroup=variable,Z_rescaled=value),paste0(nsnp_dir,"hb",bnum,"_heatmapplottedb2705neg.txt"), row.names = F, quote = F, sep = ",")

}
