# Courtney Smith - Main Analysis
# Started 3-14-2023
# Make a barplot with breakdown of trait categories for the FinnGen traits/classification, then manual classification of HLA-associated triats

library(dplyr)
library(ggplot2)

## FinnGen trait classification breakdown
an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories.txt") # Trait categories

gan <- data.table::fread("/oak/stanford/groups/pritch/users/courtrun/projects/hla/scripts/enrichment/enrich_trait_group_annotations.txt")
an <- left_join(an,gan,by=c("trait"))

count <- an %>% group_by(group) %>% count()
count <- as.data.frame(count %>% arrange(group))
count$order <- 1:nrow(count)
count$group <- ifelse(count$group=="Rheumatologic","Rheumatic",ifelse(count$group=="Rheumatologic Comorbidities","Rheumatic Comorbidities",count$group)

# Plot
p <- ggplot(count %>% filter(n>=3),aes(y=forcats::fct_reorder(group,-order),x=n,fill=group)) + #
  geom_bar(stat="identity") +
  ylab("Trait Group") +
  xlab("Number of Traits") +
  labs(fill="Trait Group")+
  xlim(0,max(count$n)+10)+
  #geom_text(aes(label=n),size=8, position=position_dodge(width=0.9), hjust=-0.25)+
  #geom_label_repel(aes(x=group,y=n,label=paste0(n),fill="white"),size=5, show.legend = FALSE) + #, min.segment.length=unit(0,'lines'),force=10,segment.size=0.25)+
  theme_classic()+
  #scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Infection","Cardiometabolic","MSK","Neoplasm","Neuro", "Organ","Other"))+
  theme(text = element_text(size = 30),legend.position="none",axis.title.y=element_blank())#+  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/trait_categories.png",width=750,height=680)
p
dev.off()

##### Plot with manual annotations
an <- data.table::fread("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories_v2.tsv") %>%
select(-Organ)# Trait categories

# Munge
an$Category <- ifelse(an$Category=="Infection","Infectious",
      ifelse(an$Category=="Neoplasm","Neoplastic",ifelse(an$Category=="Neuro","Neurologic",an$Category)))
an$Subcategory <- ifelse(an$Subcategory=="Pulmonary","Lung",
      ifelse(an$Subcategory=="Heme","Blood",an$Subcategory))

an$Category <- ifelse(an$Category=="Rheumatologic","Rheumatic",an$Category)

OUTPUT_FILE="/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/R10_trait_categories_manual.tsv"
write.table(an, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

count <- an %>% group_by(Category) %>% count()
count <- as.data.frame(count %>% arrange(Category))
count$order <- 1:nrow(count)

# Plot
p <- ggplot(count,aes(y=forcats::fct_reorder(Category,-order),x=n,fill=Category)) + #
  geom_bar(stat="identity") +
  ylab("") +
  xlab("Number of Traits") +
  labs(fill="Trait Group")+
  xlim(0,max(count$n)+10)+
  #geom_text(aes(label=n),size=8, position=position_dodge(width=0.9), hjust=-0.25)+
  #geom_label_repel(aes(x=Category,y=n,label=paste0(n),fill="white"),size=5, show.legend = FALSE) + #, min.segment.length=unit(0,'lines'),force=10,segment.size=0.25)+
  theme_classic()+
  scale_fill_brewer(palette = "Paired",breaks=c("Autoimmune","Infection","Cardiometabolic","MSK","Neoplasm","Neuro", "Organ","Other"))+
  theme(text = element_text(size = 30),legend.position="none",axis.title.y=element_blank())#+  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/trait_categories_manual.png",width=750,height=680)
p
dev.off()

## Same as above but by organ block

count <- an %>% filter(Category=="Organ") %>% group_by(Subcategory) %>% count()
count <- as.data.frame(count %>% arrange(Subcategory))
count$order <- 1:nrow(count)

# Plot
p <- ggplot(count,aes(y=forcats::fct_reorder(Subcategory,-order),x=n,fill=Subcategory)) + #
  geom_bar(stat="identity") +
  ylab("") +
  xlab("Number of Traits") +
  labs(fill="Trait Group")+
  #geom_text(aes(label=n), size=5, position=position_dodge(width=0.9), hjust=-0.25)+
  theme_classic()+
  theme(axis.text=element_text(size=30),axis.title=element_text(size=35),legend.position="none")
  #theme(text = element_text(size = 30),legend.position="none")#+  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/R10_hla_results/R10_files/trait_categories_organ.png",width=480,height=480)
p
dev.off()
