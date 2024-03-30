library(readr)
library(tibble)
library(dplyr)
library(bigmemory)
library(doParallel)
library(dplyr)
library(foreach)
library(NbClust)
library(factoextra)
library(pheatmap)

file_path <- "path/to/directory"
file_list <- list.files(path = file_path, pattern = "\\.rds$", full.names = TRUE)
all_data <- lapply(file_list, readRDS) %>% bind_rows()

age_frames <- list(`1-20` = 1:20, `21-40` = 21:40, `41-60` = 41:60, `61-80` = 61:80, `81-100` = 81:100)
#check if beginning index is included in .rds file
#NEED TO INSPECT FILE BEFORE PROGRESSING
test <- readRDS("df_GSE51057.rds")
View(test)

process_dataset <- function(file_name) {
  df <- readRDS(file_name)
  df <- df[df$age <= 100, ]
  #filter out rows that include cancer phenotypes?
  df <- df[,-c(1, 3, 4)]
  return(df)
}

compute_correlations <- function(df, age_frame) {
  df_frame <- df %>% filter(age %in% age_frame)
  sapply(df_frame[, -which(names(df_frame) == 'age'), drop = FALSE],
         function(x) cor(x, df_frame$age, use = "complete.obs"))
}

all_data <- process_dataset(all_data)
#all_data now has age then cpg columns

high_corr_cpg_sites_list <- list()
high_corr_cpg_sites <- unlist(foreach(age_frame = age_frames, .combine = c) {
  compute_correlations(all_data, age_frame)
})

filtered_df <- all_data[, c(high_corr_cpg_sites)]

#create matrix
age_range <- 0:100
num_cpg_sites <- ncol(filtered_df) - 1
B <- matrix(NA, nrow = num_cpg_sites, ncol = length(age_range))
rownames(B) <- colnames(filtered_df)[-1]
colnames(B) <- as.character(age_range)

# fill matrix with average betas
for (i in 1:num_cpg_sites) {
  for (j in age_range) {
    age_samples <- filtered_df[filtered_df[, 1] == j, i + 1]
    B[i, j] <- mean(age_samples, na.rm = TRUE)
  }
}
#matrix done!

pheatmap(B,
         scale = "row",  # Scale the rows (CpG sites)
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,
         show_colnames = TRUE,
         angle_col=45,
         fontsize_col=5)

#clustering analysis not using matrix: silhouette, ward.D2 and hierarchical clustering
#PCA?
cluster_df <- filtered_df[, -c(1)]
cluster_df <- scale(cluster_df)
fviz_nbclust(cluster_df, hcut, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")


dist_rows <- dist(B, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")

x <- #number of clusters from elbow method
clusters <- cutree(hc_rows, k = x)

pheatmap(B,
         scale = "row",  # scale the rows (CpG sites)
         cluster_rows = hc_rows,  
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,
         show_colnames = TRUE,
         angle_col = 45,
         fontsize_col = 5,
         cutree_rows = x  # Specify the number of clusters for the rows
)
