# Load libraries
library(vegan)
library(dplyr)
library(tidyr)
library(readr)

# 1. Load data
sex <- read_csv("bee_sex_plant_visitation.csv")
network <- read_csv("all_bees_network_batch2_3.csv")
scent <- read_csv("average_compound_profile_per_plant_june11.csv")

# 2. Clean and prepare scent data
# Make plant names consistent
scent$Plant_species <- gsub(" ", "_", scent$Plant_species)
sex$plant_species <- gsub(" ", "_", sex$plant_species)

# 3. Merge visit data with scent profile
sex$plant_species <- gsub(" ", "_", sex$plant_species)
merged_data <- left_join(sex, scent, by = c("plant_species" = "Plant_species"))

head(merged_data)

write.csv(merged_data, "merged_data_jun14.csv")

# 4. Average scent profile per bee species and sex
bee_profiles <- merged_data %>%
  group_by(bee_species, bee_sex) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

head(bee_profiles)


write.csv(bee_profiles, "bee_profiles_jun14.csv")

# 5. Normalize each profile to sum to 1
scent_only <- bee_profiles %>% select(-bee_species, -bee_sex)
scent_rel <- t(apply(scent_only, 1, function(x) x / sum(x)))
rownames(scent_rel) <- paste(bee_profiles$bee_species, bee_profiles$bee_sex, sep = "_")

# 6. Create metadata for PERMANOVA
metadata <- bee_profiles %>%
  mutate(ID = paste(bee_species, bee_sex, sep = "_")) %>%
  select(ID, bee_species, bee_sex)

# 7. Run PERMANOVA
dist_matrix <- vegdist(scent_rel, method = "bray")
adonis_result <- adonis2(dist_matrix ~ bee_species + bee_sex, data = metadata)

# 8. Print result
print(adonis_result)


#Test bee species effect only

adonis2(dist_matrix ~ bee_species, data = metadata)

#Test bee sex effect only

adonis2(dist_matrix ~ bee_sex, data = metadata)


#Test with interaction- whether species and sex interact (i.e., the effect of sex differs by species)

adonis2(dist_matrix ~ bee_species * bee_sex, data = metadata)

#Combine sex and species manually into a new variable (if interaction is biologically meaningful)

metadata$group <- paste(metadata$bee_species, metadata$bee_sex, sep = "_")
adonis2(dist_matrix ~ group, data = metadata)

#Pairwise comparisons-----------------------------

library(pairwiseAdonis)

# A. Species
pairwise_species <- pairwise.adonis2(
  dist_matrix ~ bee_species,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)
print(pairwise_species)

# B. Sex
pairwise_sex <- pairwise.adonis2(
  dist_matrix ~ bee_sex,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)
print(pairwise_sex)

# C. Species × Sex
metadata$group <- interaction(metadata$bee_species, metadata$bee_sex)

pairwise_combo <- pairwise.adonis2(
  dist_matrix ~ group,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)
print(pairwise_combo)




#option of permanova without the average compound profiles of plants


# Load required libraries
library(tidyverse)
library(vegan)

# Step 1: Load data
data <- read_csv("bee_sex_compound_metadata.csv")

# Step 2: Filter out rows with missing species or sex
data <- data %>% filter(!is.na(Bumblebee_species), !is.na(bee_sex))

# Step 3: Create group column for uniqueness
data <- data %>%
  mutate(group = paste(gsub(" ", "_", Bumblebee_species), bee_sex, sep = "_"))

# Step 4: Pivot to create group × compound matrix
compound_matrix <- data %>%
  select(group, Compound, Abundance) %>%
  group_by(group, Compound) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Compound, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("group")

# Step 5: Normalize each row (relative preferences)
compound_matrix <- sweep(compound_matrix, 1, rowSums(compound_matrix), FUN = "/")

# Step 6: Prepare metadata
metadata <- data %>%
  select(group, Bumblebee_species, bee_sex) %>%
  distinct()

# Step 7: Compute Bray-Curtis distance matrix
dist_matrix <- vegdist(compound_matrix, method = "bray")

# Step 8: Run PERMANOVA
adonis_result <- adonis2(dist_matrix ~ Bumblebee_species + bee_sex, data = metadata)
print(adonis_result)


#Test bee species effect only

adonis2(dist_matrix ~ Bumblebee_species, data = metadata)

#Test bee sex effect only

adonis2(dist_matrix ~ bee_sex, data = metadata)


#Test with interaction- whether species and sex interact (i.e., the effect of sex differs by species)

adonis2(dist_matrix ~ Bumblebee_species * bee_sex, data = metadata)


# Post hoc: Pairwise comparisons by species
pairwise_species <- pairwise.adonis2(
  dist_matrix ~ Bumblebee_species,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"  # Use "fdr" for less conservative correction
)

print(pairwise_species)

# Post hoc: Pairwise comparisons by sex
pairwise_sex <- pairwise.adonis2(
  dist_matrix ~ bee_sex,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)

print(pairwise_sex)

metadata$combo <- interaction(metadata$Bumblebee_species, metadata$bee_sex)

pairwise_combo <- pairwise.adonis2(
  dist_matrix ~ combo,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)

print(pairwise_combo)







# Load required libraries
library(tidyverse)
library(vegan)

# Step 1: Load data
data <- read_csv("bee_sex_compound_metadata_jun15.csv")

# Step 2: Filter out rows with missing species or sex
data <- data %>% filter(!is.na(Bumblebee_species), !is.na(bee_sex))

# Step 3: Create group column for uniqueness
data <- data %>%
  mutate(group = paste(gsub(" ", "_", Bumblebee_species), Individual, sep = "_"))

# Step 4: Pivot to create group × compound matrix
compound_matrix <- data %>%
  select(group, Compound, Abundance) %>%
  group_by(group, Compound) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Compound, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("group")

# Step 5: Normalize each row (relative preferences)
compound_matrix <- sweep(compound_matrix, 1, rowSums(compound_matrix), FUN = "/")

# Step 6: Prepare metadata
metadata <- data %>%
  select(group, Bumblebee_species, bee_sex) %>%
  distinct()

# Step 7: Compute Bray-Curtis distance matrix
dist_matrix <- vegdist(compound_matrix, method = "bray")

# Step 8: Run PERMANOVA
adonis_result <- adonis2(dist_matrix ~ Bumblebee_species + bee_sex, data = metadata)
print(adonis_result)


#Test bee species effect only

adonis2(dist_matrix ~ Bumblebee_species, data = metadata)

#Test bee sex effect only

adonis2(dist_matrix ~ bee_sex, data = metadata)


#Test with interaction- whether species and sex interact (i.e., the effect of sex differs by species)

adonis2(dist_matrix ~ Bumblebee_species * bee_sex, data = metadata)

# Post hoc: Pairwise comparisons by species
pairwise_species <- pairwise.adonis2(
  dist_matrix ~ Bumblebee_species,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"  # Use "fdr" for less conservative correction
)

print(pairwise_species)

# Post hoc: Pairwise comparisons by sex
pairwise_sex <- pairwise.adonis2(
  dist_matrix ~ bee_sex,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)

print(pairwise_sex)

metadata$combo <- interaction(metadata$Bumblebee_species, metadata$bee_sex)

pairwise_combo <- pairwise.adonis2(
  dist_matrix ~ combo,
  data = metadata,
  permutations = 999,
  p.adjust.m = "bonferroni"
)

print(pairwise_combo)









# #pairwise PERMANOVA by species
# 
# library(vegan)
# library(dplyr)
# 
# # Loop through each species
# unique_species <- unique(metadata$bee_species)
# 
# # Store results
# pairwise_results <- list()
# 
# for (species in unique_species) {
#   
#   # Subset metadata and distance matrix for this species
#   sub_meta <- metadata %>% filter(bee_species == species)
#   
#   # Skip species with only one sex represented
#   if (length(unique(sub_meta$bee_sex)) < 2) {
#     next
#   }
#   
#   # Subset the distance matrix to just the rows/columns for this species
#   sub_ids <- sub_meta$ID
#   sub_dist <- as.dist(as.matrix(dist_matrix)[sub_ids, sub_ids])
#   
#   # Run PERMANOVA within species (testing effect of sex)
#   result <- adonis2(sub_dist ~ bee_sex, data = sub_meta)
#   
#   # Store result with species name
#   pairwise_results[[species]] <- result
# }
# 
# # Print results
# for (species in names(pairwise_results)) {
#   cat("\nSpecies:", species, "\n")
#   print(pairwise_results[[species]])
# }
# 
# #Script for NMDS by Species & Sex
# 
# # Load necessary libraries
# library(vegan)
# library(ggplot2)
# library(dplyr)
# 
# # Run NMDS on full scent matrix
# set.seed(123)  # for reproducibility
# nmds <- metaMDS(scent_rel, distance = "bray", k = 2, trymax = 100)
# 
# # Extract NMDS coordinates
# nmds_points <- as.data.frame(nmds$points)
# nmds_points$ID <- rownames(nmds_points)
# 
# # Merge with metadata
# nmds_data <- left_join(nmds_points, metadata, by = "ID")
# 
# # Plotting function per species
# plot_nmds_species <- function(species_name, data = nmds_data) {
#   species_data <- filter(data, bee_species == species_name)
#   
#   ggplot(species_data, aes(x = MDS1, y = MDS2, color = bee_sex)) +
#     geom_point(size = 4, alpha = 0.8) +
#     stat_ellipse(level = 0.68, linetype = "dashed", type = "norm") +
#     labs(
#       title = paste("NMDS of Floral Scent Preferences –", species_name),
#       color = "Sex"
#     ) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(face = "bold", size = 14),
#       legend.position = "right"
#     )
# }
# 
# # Example: plot NMDS for one species
# plot_nmds_species("Bombus pratorum")
# 
# # Optional: Plot for all species in a grid (if you want multiple)
# library(patchwork)
# unique_species <- unique(nmds_data$bee_species)
# plots <- lapply(unique_species, plot_nmds_species)
# wrap_plots(plots)
# 
# 
# 
# #OPTION 2---JUN-14-------------------------
# 
# # Load required libraries
# library(vegan)
# library(dplyr)
# library(ggplot2)
# 
# # 1. Load your data
# data <- read.csv("merged_data_jun14.csv")
# 
# # 2. Prepare metadata
# metadata <- data %>%
#   select(bee_species, bee_sex) %>%
#   mutate(ID = paste0("row", row_number()))
# 
# # 3. Extract scent profile columns (assuming all numeric scent data from column 5 onwards)
# scent_data <- data[, 5:ncol(data)]
# rownames(scent_data) <- metadata$ID
# 
# scent_data[is.na(scent_data)] <- 0
# 
# scent_data <- na.omit(scent_data)
# metadata <- metadata[rownames(metadata) %in% rownames(scent_data), ]
# 
# 
# # 4. Calculate Bray–Curtis distance
# dist_matrix <- vegdist(scent_data, method = "bray")
# 
# # 5. Add ID column to metadata
# metadata$bee_species <- as.factor(metadata$bee_species)
# metadata$bee_sex <- as.factor(metadata$bee_sex)
# rownames(metadata) <- metadata$ID
# 
# # 6. PERMANOVA: interaction between species and sex
# adonis_interaction <- adonis2(dist_matrix ~ bee_species * bee_sex, data = metadata)
# print("== PERMANOVA: bee_species * bee_sex ==")
# print(adonis_interaction)
# 
# # 7. Pairwise PERMANOVA: M vs. F within each species
# unique_species <- unique(metadata$bee_species)
# pairwise_results <- list()
# 
# for (sp in unique_species) {
#   subset_meta <- metadata %>% filter(bee_species == sp)
#   if (length(unique(subset_meta$bee_sex)) < 2) next  # skip if only one sex present
#   subset_dist <- as.dist(as.matrix(dist_matrix)[subset_meta$ID, subset_meta$ID])
#   res <- adonis2(subset_dist ~ bee_sex, data = subset_meta)
#   pairwise_results[[sp]] <- res
# }
# 
# # Print pairwise results
# cat("\n== Pairwise PERMANOVA: M vs. F within species ==\n")
# for (sp in names(pairwise_results)) {
#   cat("\nSpecies:", sp, "\n")
#   print(pairwise_results[[sp]])
# }
# 
# # 8. NMDS visualization
# set.seed(42)
# nmds <- metaMDS(scent_data, distance = "bray", k = 2, trymax = 100)
# nmds_df <- as.data.frame(nmds$points)
# nmds_df$ID <- rownames(nmds_df)
# 
# # Join NMDS with metadata
# plot_data <- left_join(nmds_df, metadata, by = "ID")
# 
# # Plot NMDS
# ggplot(plot_data, aes(x = MDS1, y = MDS2, color = bee_sex, shape = bee_species)) +
#   geom_point(size = 3, alpha = 0.8) +
#   stat_ellipse(aes(group = interaction(bee_species, bee_sex)), linetype = "dashed", alpha = 0.4) +
#   labs(title = "NMDS of Bumblebee Floral Scent Preferences",
#        color = "Sex", shape = "Species") +
#   theme_minimal()
# 
# #option 3--jun14------------------------
# 
# # Load libraries
# library(vegan)
# library(dplyr)
# library(tidyr)
# library(readr)
# 
# # 1. Load data
# sex <- read_csv("bee_sex_plant_visitation.csv")
# network <- read_csv("all_bees_network_batch2_3.csv")
# scent <- read_csv("average_compound_profile_per_plant_june11.csv")
# 
# # 2. Clean plant species names to match across datasets
# scent$Plant_species <- gsub(" ", "_", scent$Plant_species)
# sex$plant_species <- gsub(" ", "_", sex$plant_species)
# 
# # 3. Merge visit data with scent profile
# merged_data <- left_join(sex, scent, by = c("plant_species" = "Plant_species"))
# 
# # 4. Count visits per bee_species, bee_sex, and plant_species
# visit_counts <- merged_data %>%
#   group_by(bee_species, bee_sex, plant_species) %>%
#   summarise(n_visits = n(), .groups = "drop")
# 
# # 5. Join visit counts with compound data
# compound_cols <- names(merged_data)[!(names(merged_data) %in% c("bee_species", "bee_sex", "plant_species"))]
# 
# merged <- visit_counts %>%
#   left_join(scent, by = c("plant_species" = "Plant_species"))
# 
# # 6. Ensure compound columns are numeric
# merged[compound_cols] <- lapply(merged[compound_cols], function(x) as.numeric(as.character(x)))
# 
# # 7. Multiply each compound by number of visits
# merged[compound_cols] <- merged[compound_cols] * merged$n_visits
# 
# # 8. Sum compound totals per bee_species × bee_sex
# bee_profiles_weighted <- merged %>%
#   group_by(bee_species, bee_sex) %>%
#   summarise(across(all_of(compound_cols), sum, na.rm = TRUE), .groups = "drop")
# 
# # Optional: Normalize rows to relative scent preference
# bee_profiles_normalized <- bee_profiles_weighted
# bee_profiles_normalized[compound_cols] <- sweep(
#   bee_profiles_normalized[compound_cols],
#   1,
#   rowSums(bee_profiles_normalized[compound_cols]),
#   FUN = "/"
# )
# 
# # Save output if needed
# write.csv(bee_profiles_normalized, "bee_profiles_relative_weighted_june14.csv", row.names = FALSE)
# 
# # 9. Prepare for PERMANOVA
# 
# # A. Build metadata dataframe
# metadata <- bee_profiles_normalized %>%
#   select(bee_species, bee_sex)
# 
# # B. Extract only compound data (as matrix)
# compound_matrix <- bee_profiles_normalized %>%
#   select(all_of(compound_cols)) %>%
#   as.matrix()
# 
# # C. Compute Bray-Curtis distance
# dist_matrix <- vegdist(compound_matrix, method = "bray")
# 
# # 10. Run PERMANOVA
# 
# # A. Test for bee species effect
# adonis_species <- adonis2(dist_matrix ~ bee_species, data = metadata)
# print(adonis_species)
# 
# # B. Test for bee sex effect
# adonis_sex <- adonis2(dist_matrix ~ bee_sex, data = metadata)
# print(adonis_sex)
# 
# # Optional: Test interaction (species × sex)
# adonis_interaction <- adonis2(dist_matrix ~ bee_species * bee_sex, data = metadata)
# print(adonis_interaction)
# 
