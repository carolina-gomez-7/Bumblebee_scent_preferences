# Load necessary libraries
library(tidyverse)
library(vegan)
library(reshape2)

# Load data
interaction_matrix <- read.csv("all_bees_network_batch2_jun21.csv", row.names = 1)
bb_compound_matrix <- read.csv("bee_compound_matrix_june11.csv", row.names = 1)
plant_compound_matrix <- read.csv("average_compound_profile_per_plant_june25.csv", row.names = 1)

# Transpose interaction matrix to get bumblebees as rows and plants as columns
interaction_matrix <- t(as.matrix(interaction_matrix))  # Bumblebee x Plant

# Step 1: Calculate total interactions
Ai <- rowSums(t(interaction_matrix))  # bumblebee totals
Aj <- colSums(interaction_matrix)     # plant totals
m <- sum(interaction_matrix)

# Step 2: Calculate expected interaction (neutral model)
expected_matrix <- outer(Ai, Aj, FUN = function(x, y) x * y / m)

# Step 3: Calculate link temperatures
link_temp_matrix <- sweep((interaction_matrix - expected_matrix), 1, Ai, "/")  # Divide row-wise by Ai

# Step 3: Calculate link temperatures
observed_matrix <- as.matrix(interaction_matrix)
link_temp_matrix <- (observed_matrix - expected_matrix) / Ai

# Step 4: Create similarity matrix between bumblebee and plants via scent
# Multiply scent matrices to get a score of 'preference match'
preference_score <- as.matrix(bb_compound_matrix) %*% t(as.matrix(plant_compound_matrix))
preference_score <- preference_score / max(preference_score)

# Optional: Correlate link temperature and preference score
temp_vector <- as.vector(link_temp_matrix)
score_vector <- as.vector(preference_score)
cor.test(temp_vector, score_vector)

# Step 5: Visualization
heatmap(link_temp_matrix, Rowv = NA, Colv = NA, scale = "none", 
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Link Temperatures (Bumblebee vs Plant)")





#option 2

# Load libraries
library(tidyverse)
library(vegan)
library(reshape2)

# Load and process matrices
# 1. Plant × Bumblebee matrix → transpose to Bumblebee × Plant
interaction_raw <- read.csv("all_bees_network_batch2_jun21.csv", row.names = 1)
interaction_matrix <- t(as.matrix(interaction_raw))  # Bumblebee x Plant

# 2. Bumblebee × Compound matrix (should be 7 x n compounds)
bb_compound_matrix <- as.matrix(read.csv("bee_compound_matrix_june11.csv", row.names = 1))

# 3. Plant × Compound matrix (should be 28 x n compounds)
plant_compound_matrix <- as.matrix(read.csv("average_compound_profile_per_plant_june25.csv", row.names = 1))

# Confirm dimensions
print(dim(interaction_matrix))       # Expect 7 x 28
print(dim(bb_compound_matrix))       # Expect 7 x n
print(dim(plant_compound_matrix))    # Expect 28 x n

# Step 1: Calculate total visits (for expected matrix)
Ai <- rowSums(interaction_matrix)  # bumblebee totals (length 7)
Aj <- colSums(interaction_matrix)  # plant totals (length 28)
m <- sum(interaction_matrix)       # total visits

# Step 2: Create expected interaction matrix
expected_matrix <- outer(Ai, Aj, FUN = function(x, y) x * y / m)  # 7 x 28

# Step 3: Calculate link temperature matrix
# Ensure Ai is a column vector for sweep
link_temp_matrix <- sweep((interaction_matrix - expected_matrix), 1, Ai, FUN = "/")  # 7 x 28

# Step 4: Preference score: Bumblebee x Plant = Bumblebee x Compounds %*% t(Plant x Compounds)
preference_score <- bb_compound_matrix %*% t(plant_compound_matrix)
preference_score <- preference_score / max(preference_score)  # Normalize to 0–1

# Optional: Correlation
cor_test <- cor.test(as.vector(link_temp_matrix), as.vector(preference_score))
print(cor_test)

# Step 5: Visualize link temperature
heatmap(link_temp_matrix, Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Link Temperatures (Bumblebee × Plant)")


#Network Visualization Script

# Load necessary libraries
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

# Prepare data for network
link_df <- as.data.frame(as.table(link_temp_matrix))  # Converts matrix to long format
colnames(link_df) <- c("Bumblebee", "Plant", "LinkTemp")

# Optional: filter out low interactions or NA
link_df <- link_df %>% filter(!is.na(LinkTemp) & LinkTemp != 0)

# Create graph object
graph <- graph_from_data_frame(link_df, directed = FALSE)

# Convert to tidygraph for ggraph plotting
graph_tidy <- as_tbl_graph(graph)

# Plot
ggraph(graph_tidy, layout = "fr") +  # 'fr' = force-directed layout
  geom_edge_link(aes(edge_alpha = abs(LinkTemp),
                     edge_width = abs(LinkTemp),
                     edge_colour = LinkTemp),
                 show.legend = TRUE) +
  geom_node_point(size = 5, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Link Temperatures in Bumblebee–Plant Network",
       edge_colour = "Link Temp") +
  theme_void()


#Script for Bipartite Network Layout

# Load libraries
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

# Step 1: Prepare data
link_df <- as.data.frame(as.table(link_temp_matrix))
colnames(link_df) <- c("Bumblebee", "Plant", "LinkTemp")

# Filter non-zero interactions
link_df <- link_df %>% filter(!is.na(LinkTemp) & LinkTemp != 0)

# Step 2: Create bipartite node list
bee_nodes <- unique(link_df$Bumblebee)
plant_nodes <- unique(link_df$Plant)

nodes <- data.frame(name = c(bee_nodes, plant_nodes),
                    type = c(rep(FALSE, length(bee_nodes)), rep(TRUE, length(plant_nodes))))  # TRUE = bumblebee

# Step 3: Create graph object
graph <- graph_from_data_frame(d = link_df, vertices = nodes, directed = FALSE)
graph_tidy <- as_tbl_graph(graph)

# Step 4: Plot in bipartite layout
ggraph(graph_tidy, layout = "bipartite") +
  geom_edge_link(aes(edge_colour = LinkTemp, edge_width = abs(LinkTemp), edge_alpha = abs(LinkTemp)),
                 show.legend = TRUE) +
  geom_node_point(aes(color = type), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, vjust = 1.5) +
  scale_edge_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_color_manual(values = c("darkgreen", "purple"), labels = c("Plant", "Bumblebee")) +
  labs(title = "Bipartite Network of Link Temperatures",
       edge_colour = "Link Temp") +
  theme_void()


#link temperature network customized

# Load libraries
library(tidyverse)
library(tibble)


# Load interaction matrix (plants × bumblebees) and transpose
interaction_raw <- read_csv("all_bees_network_batch2_jun21.csv")
interaction_matrix <- t(as.matrix(interaction_raw[,-1]))
rownames(interaction_matrix) <- colnames(interaction_raw)[-1]
colnames(interaction_matrix) <- interaction_raw[[1]]

write.csv(link_temp_matrix, "link_temp_matrix.csv")

# Read in the CSV file
link_temp_df <- read_csv("link_temp_matrix.csv")

# Convert the first column to row names
link_temp_matrix <- link_temp_df %>%
  column_to_rownames(var = names(link_temp_df)[1]) %>%
  as.matrix()


# Load link temperature matrix
link_temp_matrix <- read_csv("link_temp_matrix.csv", row.names = 1) %>% as.matrix()
link_temp_long <- as.data.frame(as.table(link_temp_matrix))
colnames(link_temp_long) <- c("Species", "Plant_species", "LinkTemp")

# Define species and plant order
plant_species_order <- c("Ranunculus acris","Tilia sp2","Epilobium angustifolium","Veronica spicata",
                         "Rhinanthus minor","Nepeta racemosa","Lavandula angustifolia","Stachys palustris",
                         "Ballota nigra","Buddleja davidii", "Ligustrum vulgare","Symphoricarpos albus",
                         "Hydrangea paniculata", "Heracleum sphondylium","Jacobaea vulgaris","Centaurea nigra",
                         "Cirsium vulgare","Cirsium arvense","Arctium minus","Rubus fruticosus agg.",
                         "Rubus caesius", "Trifolium repens","Lotus corniculatus", "Hypericum perforatum",
                         "Geranium pratense","Impatiens glandulifera", "Cornus sanguinea", "Calystegia silvatica")

species_order <- c("B_pascuorum", "B_terrestris", "B_lapidarius", "B_pratorum", "B_lucorum", "B_hortorum", "B_hypnorum")

# Prepare data
bmat <- read_csv("all_bees_network_batch2_jun21.csv")
blong <- bmat %>%
  pivot_longer(starts_with("B_"), names_to = "Species", values_to = "Weight") %>%
  group_by(Plant_species) %>%
  mutate(num_obs = sum(Weight), num_spec = nlevels(as.factor(Species[Weight != 0]))) %>%
  arrange(desc(num_spec), desc(num_obs)) %>%
  ungroup() %>%
  mutate(Plant_species = factor(Plant_species, levels = plant_species_order),
         Species = factor(Species, levels = species_order)) %>%
  mutate(across(where(is.factor), as.numeric, .names = "l_{.col}"))

# Set up coordinates
blong_selection <- blong %>%
  filter(Weight != 0) %>%
  droplevels() %>%
  mutate(Plant_species = fct_relevel(Plant_species, plant_species_order),
         Species = fct_relevel(Species, species_order)) %>%
  mutate(across(where(is.factor), as.numeric, .names = "l_{.col}"),
         x_plant = (l_Plant_species - 1) * (max(l_Species) - 1) / (max(l_Plant_species) - 1) * 3,
         x_spec = (l_Species - 1) * 3)

mid_labels <- (max(blong_selection$l_Plant_species) - 1) * max(blong_selection$l_Species) / max(blong_selection$l_Plant_species) / 2

plant_labels <- blong_selection %>%
  group_by(Plant_species, l_Plant_species) %>%
  summarise(x = min(x_plant), y = -.1)
spec_labels <- blong_selection %>%
  group_by(Species) %>%
  summarise(x = min(x_spec), y = 1.1)

# Define colors
plant_colors <- setNames(rep("chartreuse4", length(unique(plant_labels$Plant_species))), unique(plant_labels$Plant_species))
species_colors <- setNames(c("#f050ae", "#ff595e", "#ffd900", "#57c78b", "#71d6ca", "#4d9de0", "#7776bc"), unique(spec_labels$Species))
combined_colors <- c(plant_colors, species_colors)

# Merge for plotting
blong_selection_temp <- blong_selection %>%
  select(Species, Plant_species, x_plant, x_spec, l_Plant_species, l_Species) %>%
  left_join(link_temp_long, by = c("Species", "Plant_species")) %>%
  filter(!is.na(LinkTemp))

# Plot the bipartite temperature network
p1_temp <- ggplot(blong_selection_temp) +
  geom_segment(aes(x = x_plant, y = 0, xend = x_spec, yend = 0.8,
                   size = abs(LinkTemp), color = LinkTemp), lineend = "round") +
  # geom_tile(aes(x = x_plant, y = -0.052, width = 0.06, height = 0.1, fill = Plant_species)) +
  # geom_tile(aes(x = x_spec, y = 0.84, width = 0.06, height = 0.1, fill = Species)) +
  scale_size_continuous(range = c(0.25, 3)) +
  scale_color_gradient2(low = "#280659", mid = "white", high = "#f72634", midpoint = 0) +
  scale_fill_manual(values = combined_colors) +
  theme_void(base_size = 15) +
  theme(legend.position = "right") +
  annotate("text", x = mid_labels, y = -0.7, label = "Plants", fontface = 2, size = 7) +
  annotate("text", x = mid_labels, y = 1.3, label = "Bombus Species", fontface = 2, size = 7) +
  annotate("text", x = spec_labels$x, y = spec_labels$y - .05, label = spec_labels$Species, fontface = 3) +
  annotate("text", x = plant_labels$x, y = plant_labels$y - .009, label = plant_labels$Plant_species, size = 3.5, angle = 90, hjust = 1, fontface = 3)

print(p1_temp)
