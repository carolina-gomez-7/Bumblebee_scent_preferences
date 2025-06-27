#VISITED_VS_NOT_VISITED_RANDOM_FOREST
# Load required libraries
library(tidyverse)
library(randomForest)

# -------------------------------
# 1. Load the datasets
# -------------------------------
compounds <- read.csv("average_compound_profile_per_plant_june11.csv")
visitation <- read.csv("filtered_data_8april.csv")

# -------------------------------
# 2. Choose your bumblebee species
# -------------------------------
target_bumblebee <- "Bombus pascuorum"

# -------------------------------
# 3. Identify visited / not visited plants
# -------------------------------
visited <- visitation %>%
  filter(Bumblebee_species == target_bumblebee) %>%
  pull(Plant_species) %>%
  unique()

all_plants <- unique(compounds$Plant_species)
not_visited <- setdiff(all_plants, visited)

# -------------------------------
# 4. Label the data
# -------------------------------
visited_df <- tibble(Plant_species = visited, Visited = "Yes")
not_visited_df <- tibble(Plant_species = not_visited, Visited = "No")
group_labels <- bind_rows(visited_df, not_visited_df)

# Join with compound data
compounds_labeled <- compounds %>%
  inner_join(group_labels, by = "Plant_species")

# -------------------------------
# 5. Prepare data for Random Forest
# -------------------------------
# Remove species with all-zero compound profiles
compound_matrix <- compounds_labeled %>%
  select(-Plant_species, -Visited)

valid_rows <- rowSums(compound_matrix) > 0
compound_matrix <- compound_matrix[valid_rows, ]
group_info <- compounds_labeled$Visited[valid_rows]

# Build final data frame for classification
rf_data <- compound_matrix
rf_data$Visited <- as.factor(group_info)

# -------------------------------
# 6. Train Random Forest model
# -------------------------------
set.seed(123)  # for reproducibility

rf_model <- randomForest(Visited ~ ., data = rf_data, importance = TRUE, ntree = 1000)

# Print model summary
print(rf_model)

# -------------------------------
# 7. View variable importance
# -------------------------------
importance_df <- as.data.frame(importance(rf_model))
importance_df$Compound <- rownames(importance_df)

# Sort by MeanDecreaseGini
top_compounds <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice(1:15)  # Top 15 compounds

# -------------------------------
# 8. Plot top contributing compounds
# -------------------------------
library(ggplot2)

ggplot(top_compounds, aes(x = reorder(Compound, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = paste("Top Compounds Predicting Visitation by", target_bumblebee),
    x = "Compound",
    y = "Importance (MeanDecreaseGini)"
  )


anosim_result <- anosim(compound_matrix, group_info, distance = "bray")
print(anosim_result)

library(vegan)

# Bray-Curtis dissimilarity
dist_matrix <- vegdist(compound_matrix, method = "bray")

# MRPP
mrpp_result <- mrpp(dist_matrix, group_info, distance = "bray")
print(mrpp_result)
