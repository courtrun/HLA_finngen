library("dplyr")
library("dendextend")
library("RColorBrewer")
library(ggplot2)
library(ComplexHeatmap)
library(cowplot)


# Set the dimensions of the matrix
num_rows <- 30
num_cols <- 30

# Generate the random matrix
set.seed(30)
m1 <- data.frame(matrix(sample(c(0, 1), num_rows * num_cols, replace = TRUE), nrow = num_rows, ncol = num_cols))
set.seed(50)
m2 <- data.frame(matrix(sample(c(0, 1), num_rows * num_cols, replace = TRUE), nrow = num_rows, ncol = num_cols))
m <- distinct(bind_rows(m1,m2))
m$ID <- c(paste0("ID",c(1:nrow(m))))
full <- m
m <- full %>% select(-ID)
rownames(m) <- full$ID

dis <- dist(m, method = "binary") # Compute the dissimilarity matrix
hc <- hclust(dis, method = "complete") # Perform hierarchical clustering
dend <- as.dendrogram(hc)
plot(hc)
ngroups=10
rect.hclust(hc,k=ngroups)
hc$order
cluster_labels <- cutree(hc, k = ngroups)

# Add info
m$ID <- full$ID
m$cluster <- cluster_labels
m <- left_join(m,data.frame(cluster=c(1:10),colors=brewer.pal(10, "Paired")),by=c("cluster"))

#filter(m,ID %in% c("ID15","ID35","ID47","ID55","ID30","ID9","ID50"))

# Reorder to dendrogram and get colors
m <- m[hc$order,]
color_list <- m$colors

m$dendroname <- paste0("Cluster ",m$cluster,": ",m$ID)

dend %>% set("labels_col",color_list) %>% set("labels_cex", 0.75) %>% set("branches_k_color",color_list) %>% plot(horiz=TRUE,axes=FALSE)

#color_order <- get_leaves_branches_col(dend %>% set("labels_col",color_list) %>% set("labels_cex", 0.75) %>% set("branches_k_color",color_list)) # get list of colors in order from bottom to top of sorted dendrogram
#setdiff(color_list,color_order) # the same

#####
# Randomly subset to 2 from each group
num_rows_to_select=2
set.seed(30)
subset_df <- as.data.frame(m %>% group_by(cluster) %>% arrange(ID) %>% slice(1:num_rows_to_select))

# Reformat
subset_df$dendroname <- paste0("Cluster ",subset_df$cluster,": ",subset_df$ID)
m2 <- subset_df %>% select(-ID,-cluster,-colors,-dendroname)
rownames(m2) <- subset_df$dendroname

# Cluster
dis <- dist(m2, method = "binary") # Compute the dissimilarity matrix
hc <- hclust(dis, method = "complete") # Perform hierarchical clustering
dend <- as.dendrogram(hc)
plot(hc)
ngroups=10
rect.hclust(hc,k=ngroups)
hc$order
m2$dendroname <- subset_df$dendroname
m2 <- left_join(m2,m %>% select(dendroname,cluster,colors),by=c("dendroname"))
m2 <- m2[hc$order,]
color_list <- m2$colors

d <- dend %>% set("labels_col",color_list) %>% set("labels_cex", 1) %>% set("branches_k_color",color_list) %>% sort(type="labels")
sorted_order <- unique(unlist(partition_leaves(d)))
color_order <- get_leaves_branches_col(d) # get list of colors in order from bottom to top of sorted dendrogram

png("tempdendro.png",width=672,height=672)
par(mar = c(2,2,2,10))
plot(d,horiz=TRUE,axes=FALSE)
dev.off()

melted_df <- reshape2::melt(m2 %>% select(-cluster,-colors))
melted_df <- left_join(melted_df,m2 %>% select(cluster,dendroname,colors),by=c("dendroname"))

# make sure value column is only 1 or 0
filter(melted_df,value!=0&value!=1) # should be empty
melted_df$vc <- ifelse(melted_df$value==1,melted_df$colors,"white")
melted_df$dendroname <- factor(melted_df$dendroname, levels = sorted_order)

# Plot the heatmap using ggplot2
ggplot(melted_df, aes(x = variable, y = dendroname, fill = factor(vc))) +
  geom_tile(color="black",size=1)+
  scale_fill_identity()+labs(x="SNPs",y="Haplotypes")
length(unique(melted_df$colors))
ggsave("tempheatmap.png",width=7,height=7)
