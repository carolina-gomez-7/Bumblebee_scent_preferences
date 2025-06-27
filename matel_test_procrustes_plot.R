# Load necessary libraries
library(ggplot2)
library(vegan)
library(reshape2)
library(pheatmap)

# Load the CSV files
scent_mat <- as.matrix(read.csv("bray_curtis_distances_plants_jun25.csv", row.names = 1, check.names = FALSE))
visit_mat <- as.matrix(read.csv("visit_distance_matrix.csv", row.names = 1, check.names = FALSE))

# Check if row and column names match
if (!all(rownames(scent_mat) == colnames(scent_mat)) || !all(rownames(visit_mat) == colnames(visit_mat))) {
  stop("Row and column names of the matrices must match.")
}

# Convert to dist objects (lower triangle only)
scent_dist <- as.dist(scent_mat)
visit_dist <- as.dist(visit_mat)

# Create a data frame for Mantel scatterplot
mantel_df <- data.frame(
  ScentDissimilarity = as.vector(scent_dist),
  VisitDissimilarity = as.vector(visit_dist)
)

# Plot: Scatterplot with loess smoothing
ggplot(mantel_df, aes(x = ScentDissimilarity, y = VisitDissimilarity)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", fill = "lightblue", se = TRUE) +
  labs(
    title = "Relationship Between Floral Scent and Bumblebee Visitation Dissimilarities",
    x = "Floral Scent Dissimilarity (Bray–Curtis)",
    y = "Visitation Dissimilarity (Bray–Curtis)"
  ) +
  theme_minimal()


# Compute residuals
mantel_residuals <- visit_mat - scent_mat

# Set diagonal to NA for clarity
diag(mantel_residuals) <- NA

# Plot heatmap
pheatmap(
  mantel_residuals,
  main = "Heatmap of Mantel Residuals (Visitation - Scent)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  na_col = "grey90",
  display_numbers = FALSE
)


#Procrustes Overlay Visualization


# Load required libraries
library(vegan)
library(ggplot2)

# Step 1: Load your distance matrices
# Replace with your actual file paths
scent_mat <- as.matrix(read.csv("bray_curtis_distances_plants_jun25.csv", row.names = 1, check.names = FALSE))
visit_mat <- as.matrix(read.csv("visit_distance_matrix.csv", row.names = 1, check.names = FALSE))

# Step 2: Convert matrices to distance objects
scent_dist <- as.dist(scent_mat)
visit_dist <- as.dist(visit_mat)

# Step 3: Perform PCoA
pcoa_scent <- cmdscale(scent_dist, eig = TRUE, k = 2)
pcoa_visit <- cmdscale(visit_dist, eig = TRUE, k = 2)

# Step 4: Procrustes analysis (scent onto visitation)
proc <- procrustes(pcoa_visit$points, pcoa_scent$points, symmetric = TRUE)

# Step 5: Manually extract coordinates from Procrustes object
procrustes_df <- data.frame(
  Species = rownames(scent_mat),
  Scent_X = proc$Yrot[,1],
  Scent_Y = proc$Yrot[,2],
  Visit_X = proc$X[,1],
  Visit_Y = proc$X[,2]
)

# Step 6: Plot Procrustes results
ggplot(procrustes_df) +
  geom_segment(aes(x = Scent_X, y = Scent_Y, xend = Visit_X, yend = Visit_Y),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.5, color = "gray") +
  geom_point(aes(x = Scent_X, y = Scent_Y), color = "darkgreen", size = 3) +
  geom_point(aes(x = Visit_X, y = Visit_Y), color = "darkblue", size = 3) +
  geom_text(aes(x = Visit_X, y = Visit_Y, label = Species), hjust = 1.2, size = 3) +
  labs(
    title = "Procrustes Alignment of PCoA Configurations",
    subtitle = "Floral Scent (green) vs. Bumblebee Visitation (blue)",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2"
  ) +
  theme_minimal()

protest(procrustes_df)

#PCoA for Floral Scent Dissimilarity

# Load vegan and ape for ordination
library(vegan)
library(ape)

# PCoA (Principal Coordinates Analysis) on scent distances
pcoa_scent <- cmdscale(scent_dist, eig = TRUE, k = 2)

# Create data frame for plotting
scent_pcoa_df <- data.frame(
  Species = rownames(scent_mat),
  Axis1 = pcoa_scent$points[,1],
  Axis2 = pcoa_scent$points[,2]
)

# Plot
ggplot(scent_pcoa_df, aes(x = Axis1, y = Axis2, label = Species)) +
  geom_point(color = "darkgreen", size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  labs(
    title = "PCoA of Floral Scent Dissimilarities",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2"
  ) +
  theme_minimal()

#PCoA for Bumblebee Visitation Dissimilarity

# PCoA on visitation distances
pcoa_visit <- cmdscale(visit_dist, eig = TRUE, k = 2)

# Create data frame
visit_pcoa_df <- data.frame(
  Species = rownames(visit_mat),
  Axis1 = pcoa_visit$points[,1],
  Axis2 = pcoa_visit$points[,2]
)

# Plot
ggplot(visit_pcoa_df, aes(x = Axis1, y = Axis2, label = Species)) +
  geom_point(color = "darkblue", size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  labs(
    title = "PCoA of Bumblebee Visitation Dissimilarities",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2"
  ) +
  theme_minimal()


