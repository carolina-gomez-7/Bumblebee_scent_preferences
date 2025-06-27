#TANGLEGRAM jun18
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ape)
library(dendextend)
library(stringr)
library(phytools)
library(rcompanion)
library(vegan)

# Load the dataset
file_path <- "bee_distances_mat_jun18.csv"  # Ensure the correct path
bray_curtis_distances <- as.dist(read.csv(file_path, row.names = 1))

# Load and process the phylogenetic tree
tree <- read.tree("partitions.txt.treefile")

plot(tree)

# Clean and standardize tree tip labels
tree$tip.label <- str_replace_all(tree$tip.label, "[_.0-9]+$", "")

# Root the tree at the midpoint
tree <- midpoint.root(tree)

# Check if the tree is ultrametric, binary, and rooted
if (!is.ultrametric(tree)) {
  tree <- chronos(tree)  # Convert to ultrametric if necessary
}

if (!is.binary.tree(tree)) {
  tree <- multi2di(tree)  # Convert to binary if necessary
}

if (!is.rooted(tree)) {
  tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)  # Ensure it's rooted
}

# Perform hierarchical clustering
hc_chemical <- hclust(bray_curtis_distances, method = "average")


# Convert tree directly to a dendrogram
phylo_dendro <- as.dendrogram(tree)

# # Convert phylogenetic tree to dendrogram properly
# phylo_dendro <- as.dendrogram(hclust(cophenetic(tree)))
# 
# p1 <-cophenetic(tree)
# p2 <- hclust(p1)
  
chemical_dendro <- as.dendrogram(hc_chemical)

#plot tanglegram
tanglegram(phylo_dendro, chemical_dendro,
           main_left = "Phylogenetic Tree",
           main_right = "Scent Preference Dendrogram",
           common_subtrees_color_lines = TRUE,
           highlight_distinct_edges = TRUE,
           lab.cex = 0.8,
           fast = TRUE)

# Match labels and color branches
phylo_dendro <- phylo_dendro %>% set("branches_k_color", k = 3)  # e.g. 3 clusters
chemical_dendro <- chemical_dendro %>% set("branches_k_color", k = 3)

# Plot tanglegram with coloring
tanglegram(phylo_dendro, chemical_dendro,
           main = "Tanglegram: Phylogeny vs. Scent Preference",
           common_subtrees_color_lines = TRUE,
           highlight_distinct_edges = TRUE,
           margin_inner = 7,
           lab.cex = 0.8,
           fast = TRUE)

#Compute entanglement and cophenetic correlation

entanglement_value <- entanglement(dendlist(phylo_dendro, chemical_dendro))
print(paste("Entanglement score:", round(entanglement_value, 3)))


#Cophenetic correlation (Pearson correlation between two dendrograms)

cophenetic_correlation <- cor_cophenetic(phylo_dendro, chemical_dendro)
print(paste("Cophenetic correlation:", round(cophenetic_correlation, 3)))

# Compute Pearson correlation
dend_list <- dendlist(phylo_dendro, chemical_dendro)
pearson_corr <- cor.dendlist(dend_list, method_coef = "pearson")
print(paste("Pearson Correlation: ", round(pearson_corr, 4)))


#TANGLEGRAM OPTION 2- MANTEL TEST AND BOOTSTRAP VALUES

# Load required packages
library(ape)
library(vegan)
library(ggtree)
library(tidyverse)

# --- PART 1: MANTEL TEST ---

# Load and prepare Bray–Curtis distance matrix
bray <- read.csv("bee_distances_mat_jun18.csv", row.names = 1)
bray <- as.dist(as.matrix(bray))

# Load and process the phylogenetic tree
tree <- read.tree("partitions.txt.treefile")

# Optionally clean tip labels
tree$tip.label <- gsub("(_[0-9]+|\\.[0-9]+)$", "", tree$tip.label)  # remove suffixes

# Ensure the tree is ultrametric and binary for distance calculation
if (!is.ultrametric(tree)) tree <- chronos(tree)
if (!is.binary.tree(tree)) tree <- multi2di(tree)

# Compute pairwise phylogenetic distances
phylo_dist <- cophenetic(tree)

# Match dimensions
common <- intersect(rownames(phylo_dist), labels(bray))
phylo_dist <- phylo_dist[common, common]
bray <- as.matrix(bray)[common, common]
bray <- as.dist(bray)

# Convert phylogenetic distance to dist object
phylo_dist <- as.dist(phylo_dist)

# Run Mantel test
mantel_result <- mantel(phylo_dist, bray, method = "pearson", permutations = 9999)
print(mantel_result)

# --- PART 2: PLOT TREE WITH BOOTSTRAP VALUES ---

# Visualize tree with bootstrap values (if available)
if (!is.null(tree$node.label)) {
  ggtree(tree) +
    geom_tiplab(size = 1) +
    geom_text2(aes(subset = !isTip, label = label), size = 3, hjust = -0.2) +
    ggtitle("Phylogeny with Bootstrap Support") +
    theme_tree2()
} else {
  message("Tree has no bootstrap values in node labels.")
}
