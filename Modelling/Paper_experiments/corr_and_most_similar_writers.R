setwd("C:/Users/ltzai/Desktop/PhD/Handwritten_Loop_characters/handwriting-evidence-evaluation-fourier-based-feature-extraction/Modelling")
library(dplyr)
library(readxl)

IAM_data <- read_excel("IAM_fourier_features_dataset/DB_loop_handwriting.xlsx")
IAM_data = as.data.frame(IAM_data)

library(PerformanceAnalytics)

chart.Correlation(
  R          = IAM_data[IAM_data$writer_id=='88',2:9],
  histogram  = TRUE,            # histograms on diagonal
  method     = "pearson",       # or "spearman", "kendall"
  pch        = 19               # point symbol in scatterplots
)

IAM_data = cbind(scale(IAM_data[,1:9]),IAM_data[,10:ncol(IAM_data)])
# Step 1: Rename columns for clarity

colnames(IAM_data)[10:11] <- c("writer_id", "character")

feature_cols <- colnames(IAM_data)[1:9]

# Step 3: Mean per writer-character
writer_character_means <- IAM_data %>%
  group_by(writer_id, character) %>%
  summarise(across(all_of(feature_cols), mean, na.rm = TRUE), .groups = "drop")

# Step 4: Mean per writer across characters
writer_profiles <- writer_character_means %>%
  group_by(writer_id) %>%
  summarise(across(all_of(feature_cols), mean, na.rm = TRUE), .groups = "drop")

# Step 4: Compute distance matrix between writers
# Remove writer_id column for distance computation
feature_matrix <- writer_profiles %>%  dplyr::select(-writer_id)
writer_ids <- writer_profiles$writer_id

# Use Euclidean distance (you can change to e.g., "manhattan" or "cosine" via `proxy::dist`)
dists <- dist(feature_matrix, method = "euclidean")

D <- as.matrix(dists)
diag(D) <- NA

# work only with upper triangle to avoid duplicates
D_upper <- D
D_upper[lower.tri(D_upper, diag = TRUE)] <- NA

k <- which.min(D_upper)
pair_idx <- arrayInd(k, dim(D_upper))

closest_pair_ids <- writer_ids[pair_idx]
closest_pair_dist <- D_upper[pair_idx]


pairs_idx <- which(D_upper < 0.5, arr.ind = TRUE)

pairs_table <- data.frame(
  writer1  = writer_ids[pairs_idx[, 1]],
  writer2  = writer_ids[pairs_idx[, 2]],
  distance = D_upper[pairs_idx]
)


# Step 5: Hierarchical clustering
hc <- hclust(dists, method = "average")  # or "complete", "ward.D2", etc.

# Plot dendrogram to visually inspect
plot(hc, labels = writer_ids, main = "Writer Similarity (Hierarchical Clustering)")

# Step 6: Cut dendrogram into clusters
# Try getting clusters of ~20 or ~30 writers
# Start with number of clusters that gives groups of size ~20
cluster_labels <- cutree(hc, k = floor(length(writer_ids) / 20))  # try also 30 later

# Add cluster info to writer_profiles
writer_profiles$cluster <- cluster_labels

# Step 7: Find most cohesive cluster (lowest intra-cluster distance)
find_most_similar_k <- function(k_desired) {
  # Compute full distance matrix
  dist_mat <- as.matrix(dist(feature_matrix))
  diag(dist_mat) <- Inf  # ignore self-distances
  
  # Find all combinations of k_desired writers (too slow for large n)
  # Instead: greedy method
  
  # Start with the pair of writers with smallest distance
  min_idx <- which(dist_mat == min(dist_mat), arr.ind = TRUE)[1, ]
  selected <- c(min_idx[1], min_idx[2])
  remaining <- setdiff(1:nrow(dist_mat), selected)
  
  while (length(selected) < k_desired) {
    avg_dists <- sapply(remaining, function(i) {
      mean(dist_mat[i, selected])
    })
    next_best <- remaining[which.min(avg_dists)]
    selected <- c(selected, next_best)
    remaining <- setdiff(remaining, next_best)
  }
  
  return(writer_profiles$writer_id[selected])
}

# Get exactly 20 most similar writers
most_similar_20 <- find_most_similar_k(20)

# Get exactly 30 most similar writers
most_similar_30 <- find_most_similar_k(30)

# Output
cat("Most similar 20 writers (exactly):\n")
print(most_similar_20)

cat("Most similar 30 writers (exactly):\n")
print(most_similar_30)



dist_matrix <- as.matrix(dist(feature_matrix))  # square symmetric matrix

# Step 2: Set the diagonal to Inf so we ignore self-comparisons
diag(dist_matrix) <- Inf

# Step 3: Find indices of the smallest value (i.e., most similar pair)
min_idx <- which(dist_matrix == min(dist_matrix), arr.ind = TRUE)[1, ]

# Step 4: Get the writer IDs
most_similar_pair <- writer_profiles$writer_id[c(min_idx[1], min_idx[2])]

# Print result
cat("The two most similar writers are:", most_similar_pair[1], "and", most_similar_pair[2], "\n")

