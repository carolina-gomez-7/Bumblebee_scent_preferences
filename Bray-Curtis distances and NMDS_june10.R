# Load required packages
library(vegan)       # for vegdist and ordination
library(ggplot2)     # for clean plotting
library(tidyverse)   # for data manipulation

# Step 1: Load data
# Replace with your actual file paths
scent <- read.csv("average_compound_profile_per_plant_june11.csv", row.names = 1, check.names = FALSE)
visit <- read.csv("all_bees_network_batch2.csv", row.names = 1, check.names = FALSE)

# Step 2: Standardize plant names (replace underscores with spaces)
rownames(visit) <- gsub("_", " ", rownames(visit))

# Step 3: Keep only common plant species
common_plants <- intersect(rownames(scent), rownames(visit))
scent <- scent[common_plants, ]
visit <- visit[common_plants, ]

# Step 4: Create bee × plant matrix (transpose)
bee_plant <- t(visit)

# Step 5: Multiply bee × plant with plant × compound → bee × compound matrix
bee_compound <- as.matrix(bee_plant) %*% as.matrix(scent)

write.csv(bee_compound, "bee_compound_matrix_june11.csv")

# Step 6: Normalize rows (relative preferences)
bee_compound_rel <- sweep(bee_compound, 1, rowSums(bee_compound), FUN = "/")

# Step 7: Remove bees with NA or zero rows
bee_compound_rel <- bee_compound_rel[complete.cases(bee_compound_rel), ]
bee_compound_rel <- bee_compound_rel[rowSums(bee_compound_rel) > 0, ]

# Step 8: Compute Bray-Curtis distances between bee species
bee_dist <- vegdist(bee_compound_rel, method = "bray")

bee_dist_mat <- as.matrix(bee_dist)

write.csv(bee_dist_mat, "bee_distances_mat_jun11.csv")

# Step 9: Perform MDS (PCoA also works)
mds <- cmdscale(bee_dist, k = 2)
mds_df <- as.data.frame(mds)
mds_df$bee <- rownames(mds)

# Step 10: Plot
ggplot(mds_df, aes(x = V1, y = V2, label = bee)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(vjust = -0.7, size = 3) +
  labs(
    title = "Clustering of Bumblebee Species by Floral Scent Preferences",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  ) +
  theme_minimal()

#option plotting envfit

# Assume you have:
# - mds: your ordination result (from cmdscale or metaMDS)
# - bee_compound_rel: matrix of bee × compound relative preferences

# Run envfit
fit <- envfit(mds, bee_compound_rel, permutations = 999)

# View significant compounds
print(fit)

# Plot ordination with vectors
plot(mds, type = "n")
points(mds, pch = 19, col = "blue")
plot(fit, p.max = 0.05, col = "red")  # Only show significant vectors


#option plotting NMDS

# Load required packages
library(vegan)
library(ggplot2)

# Step 1: Prepare your bee × compound relative preference matrix
# Assuming you've already created this:
# bee_compound_rel = normalized matrix (rows = bee species, cols = compounds)

# Remove rows with zero total (if not done yet)
bee_compound_rel <- bee_compound_rel[rowSums(bee_compound_rel) > 0, ]

# Step 2: Run NMDS
set.seed(123)  # for reproducibility
nmds <- metaMDS(bee_compound_rel, distance = "bray", k = 2, trymax = 100)

# Check stress value
print(nmds$stress)

# Step 3: Fit compound vectors using envfit
fit <- envfit(nmds, bee_compound_rel, permutations = 999)

# Step 4: Base R plot with vectors
plot(nmds, type = "n", main = "NMDS of Bumblebee Floral Scent Preferences", 
     xlim = c(-0.5, 0.5),
     ylim = c(-0.5, 0.5))
points(nmds, display = "sites", pch = 21, bg = "steelblue", cex = 2)
text(nmds, display = "sites", pos = 3, cex = 0.8)
plot(fit, p.max = 0.05, col = "red")  # Only significant compound vectors

# Optional: print vector statistics
print(fit)


# Optional: Get coordinates for ggplot
site_scores <- scores(nmds, display = "sites")
vector_scores <- scores(fit, display = "vectors")
vector_pvals <- fit$vectors$pvals
significant_vectors <- vector_scores[vector_pvals <= 0.05, , drop = FALSE]


#Final NMDS + envfit script with scaled arrows and axis limits

# Load required package
library(vegan)

# Step 1: Prepare bee × compound relative preference matrix
bee_compound_rel <- bee_compound_rel[rowSums(bee_compound_rel) > 0, ]

# Step 2: Run NMDS
set.seed(123)
nmds <- metaMDS(bee_compound_rel, distance = "bray", k = 2, trymax = 100)
print(paste("NMDS stress:", round(nmds$stress, 4)))

# Step 3: Fit compound vectors using envfit
fit <- envfit(nmds, bee_compound_rel, permutations = 999)

# Step 4: Scale vectors
scaling_factor <- 0.5  # Adjust this value to control arrow length
vectors_scaled <- scores(fit, display = "vectors") * scaling_factor

# Step 5: Plot NMDS with axis limits and scaled vectors
plot(nmds,
     type = "n",
     main = "NMDS of Bumblebee Floral Scent Preferences",
     xlim = c(-0.5, 0.5),
     ylim = c(-0.5, 0.5))

# Add bee species points
points(nmds, display = "sites", pch = 21, bg = "steelblue", cex = 2.5)
text(nmds, display = "sites", pos = 3, cex = 0.8)

# Add arrows (scaled vectors)
arrows(
  0, 0,
  vectors_scaled[, 1], vectors_scaled[, 2],
  length = 0.1, col = "red", lwd = 1.5
)

# Add labels to the vector tips
text(
  vectors_scaled[, 1], vectors_scaled[, 2],
  labels = rownames(vectors_scaled),
  col = "red", cex = 0.7, pos = 4
)


#NMDS of Bumblebee Scent Preferences with Significant Vectors Only


# Load required package
library(vegan)

# Step 1: Prepare bee × compound relative preference matrix
bee_compound_rel <- bee_compound_rel[rowSums(bee_compound_rel) > 0, ]

# Step 2: Run NMDS
set.seed(123)
nmds <- metaMDS(bee_compound_rel, distance = "bray", k = 2, trymax = 100)
cat("NMDS stress:", round(nmds$stress, 4), "\n")

# Step 3: Fit compound vectors using envfit
fit <- envfit(nmds, bee_compound_rel, permutations = 999)

# Step 4: Extract significant vectors only (p ≤ 0.05)
pvals <- fit$vectors$pvals
significant_vectors <- scores(fit, display = "vectors")[pvals <= 0.05, , drop = FALSE]

# Step 5: Scale vector length
scaling_factor <- 0.5
vectors_scaled <- significant_vectors * scaling_factor

# Step 6: Plot NMDS with axis limits and only significant vectors
plot(nmds,
     type = "n",
     main = "NMDS of Bumblebee Floral Scent Preferences\n(Significant Compounds Only)",
     xlim = c(-0.2, 0.2),
     ylim = c(-0.5, 0.5))

# Add bumblebee species points
points(nmds, display = "sites", pch = 21, bg = "steelblue", cex = 3)
text(nmds, display = "sites", pos = 3, cex = 0.8)

# Plot only significant vectors
arrows(
  0, 0,
  vectors_scaled[, 1], vectors_scaled[, 2],
  length = 0.05, col = "red", lwd = 1.5
)

# Label significant compound vectors
text(
  vectors_scaled[, 1], vectors_scaled[, 2],
  labels = rownames(vectors_scaled),
  col = "red", cex = 0.6, pos = 4
)


#version with ggrepel to avoid overlaping

# Load libraries
library(vegan)
library(ggplot2)
library(ggrepel)

# Step 1: Run NMDS and envfit (same as before)
bee_compound_rel <- bee_compound_rel[rowSums(bee_compound_rel) > 0, ]
set.seed(123)
nmds <- metaMDS(bee_compound_rel, distance = "bray", k = 2, trymax = 100)
fit <- envfit(nmds, bee_compound_rel, permutations = 999)

# Step 2: Extract scores
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$bee <- rownames(site_scores)

pvals <- fit$vectors$pvals
vectors <- scores(fit, display = "vectors")
sig_vectors <- vectors[pvals <= 0.05, , drop = FALSE]

# Step 3: Prepare vector data for ggplot
scaling_factor <- 0.3
vector_df <- as.data.frame(sig_vectors * scaling_factor)
vector_df$compound <- rownames(vector_df)
colnames(vector_df)[1:2] <- c("NMDS1", "NMDS2")

# Step 4: Plot using ggplot2 + ggrepel
ggplot(site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text_repel(aes(label = bee), size = 3.5, color = "black") +
  geom_segment(data = vector_df,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "red", size = 1) +
  geom_text_repel(data = vector_df,
                  aes(x = NMDS1, y = NMDS2, label = compound),
                  color = "red", size = 3, max.overlaps = Inf, segment.color = "grey40") +
  labs(
    title = "NMDS of Bumblebee Species by Floral Scent Preferences",
    subtitle = "Only significant scent compounds (envfit, p ≤ 0.05)",
    x = "NMDS1", y = "NMDS2"
  ) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  theme_minimal()



# Step 5 (Optional): ggplot version of NMDS + vectors
library(ggrepel)

site_df <- as.data.frame(site_scores)
site_df$bee <- rownames(site_df)

vector_df <- as.data.frame(significant_vectors)
vector_df$compound <- rownames(vector_df)
# Rename columns for plotting
colnames(vector_df)[1:2] <- c("NMDS1", "NMDS2")


ggplot(site_df, aes(x = NMDS1, y = NMDS2, label = bee)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text_repel(size = 3) +
  geom_segment(data = vector_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "red") +
  geom_text_repel(data = vector_df, aes(x = NMDS1, y = NMDS2, label = compound),
                  color = "red", size = 3, segment.color = "red") +
  labs(title = "NMDS of Bumblebee Species by Floral Scent Preferences",
       x = "NMDS1", y = "NMDS2") +
  theme_minimal()


ggplot(site_df, aes(x = NMDS1, y = NMDS2, label = bee)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text_repel(size = 3) +
  geom_segment(data = vector_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, color = compound_class),
               arrow = arrow(length = unit(0.25, "cm")), size = 1) +
  geom_text_repel(data = vector_df, aes(x = NMDS1, y = NMDS2, label = compound, color = compound_class),
                  size = 3, segment.color = "grey40") +
  labs(
    title = "NMDS of Bumblebee Species by Floral Scent Preferences",
    x = "NMDS1", y = "NMDS2", color = "Compound Class"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Terpenoid" = "forestgreen", "Aromatic" = "darkred", "Other" = "gray40"))


