library(dplyr)
library(ggplot2)

melted_df1 <- data.table::fread(paste0("/home/ivm/haplo_analysis/snps1000/hb",1,"_heatmapplotted.txt"))
melted_df2 <- data.table::fread(paste0("/home/ivm/haplo_analysis/snps1000/hb",2,"_heatmapplotted.txt"))
melted_df3 <- data.table::fread(paste0("/home/ivm/haplo_analysis/snps1000/hb",3,"_heatmapplotted.txt"))

plot_zcut=3

melted_df_comb <- bind_rows(melted_df1 %>% mutate(block=1),melted_df2 %>% mutate(block=2),melted_df3 %>% mutate(block=3)) %>%
mutate(direction=ifelse(value > plot_zcut,"above2",ifelse(value < (-plot_zcut),"below2","neither")))

updown <- melted_df_comb %>% group_by(variable,block) %>%
filter(any(direction=="above2")&any(direction=="below2")) %>% filter(direction!="neither") %>%
arrange(block,variable)

data.long <- updown %>% group_by(direction,variable,block) %>% count() %>% arrange(block,variable,-n)
data.long$n <- ifelse(data.long$direction=="above2",data.long$n,-1*data.long$n)

# Plot
ggplot(data.long, aes(x = variable, y = n, fill = direction)) + facet_wrap(~block,scales="free_y")+
geom_bar(data = subset(data.long, direction == "above2"), stat="identity") +
geom_bar(data = subset(data.long, direction == "below2"), stat="identity") +
coord_flip() +theme(axis.text=element_text(size=35),
                    axis.title = element_text(size=40),
                    legend.text=element_text(size=35),
                    legend.title=element_text(size=40))+
labs(fill = "Z-Score", y = "Number of Trait Associations", x = "Haplotype Group") +
theme_minimal()+
scale_fill_discrete(labels=c(paste0("Z > ",plot_zcut),paste0("Z < -",plot_zcut)))+
scale_x_continuous(breaks = seq(0, max(data.long$variable), by = 1))
ggsave(paste0(fdir,"snps",nsnp,"/figures_askout07012024/hapupdown.tiff"),height=10,width=12)

# Explore Results
data.long %>% group_by(direction) %>% summarize(mean(abs(n)))
head(data.long)
data.long %>% group_by(direction,variable,block) %>% summarize(total_n=sum(n)) %>% arrange(block,-total_n)
updown %>% group_by(direction,variable,block) %>% count() %>% arrange(block,variable,-n)
head(updown)

filter(updown,variable=="7"&block=="1")
as.data.frame(filter(updown,variable=="7"&block=="1")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="22"&block=="1")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="14"&block=="1")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="6"&block=="1")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="22"&block=="1"))
as.data.frame(filter(updown,variable=="10"&block=="2")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="3"&block=="2")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="5"&block=="2")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="7"&block=="2")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="10"&block=="2"))
as.data.frame(filter(updown,variable=="22"&block=="3")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="20"&block=="3")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="16"&block=="3")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="15"&block=="3")) %>% group_by(direction) %>% count()
as.data.frame(filter(updown,variable=="14"&block=="3")) %>% group_by(direction) %>% count()

# Exploring all in one direction
melted_df_comb %>% filter(variable==49 & abs(value) > plot_zcut)
melted_df_comb %>% filter(variable==49 & value > plot_zcut)
