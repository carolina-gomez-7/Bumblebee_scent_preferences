# Load libraries
library(tidyverse)
library(vegan)
library(reshape2)

# Load the data
df <- read.csv("bee_sex_compound_normalised_matrix_2.csv", row.names = 1)

# Extract bee species and sex from row names
df <- df %>%
  rownames_to_column("bee_id") %>%
  separate(bee_id, into = c("species", "sex"), sep = "_(?=[MF]$)", remove = FALSE)

colnames(df)


# Pivot data to wide form for paired Wilcoxon tests
paired_list <- df %>%
  pivot_longer(cols = -c(bee_id, species, sex), names_to = "compound", values_to = "abundance") %>%
  pivot_wider(names_from = sex, values_from = abundance) %>%
  filter(!is.na(F) & !is.na(M))

paired_list <- read.csv("paired_list_jun15.csv")

# Perform Wilcoxon signed-rank test per species
results <- paired_list %>%
  group_by(bee_species) %>%
  summarise(
    p_value = tryCatch({
      wilcox.test(F, M, paired = TRUE)$p.value
    }, error = function(e) NA_real_)
  )

# Print p-values
print(results)

# Optional: Adjust p-values for multiple testing (FDR)
results$p_adjusted <- p.adjust(results$p_value, method = "fdr")

print(results)

# --- Optional visualization: boxplot of differences ---
paired_diff <- paired_list %>%
  mutate(diff = F - M)

ggplot(paired_diff, aes(x = bee_species, y = diff)) +
  geom_boxplot(aes(fill = bee_species)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Sex Differences in Floral Scent Preferences",
    y = "Female - Male Relative Abundance",
    x = "Bumblebee Species"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip()


#Barplot of Adjusted P-values per Species
ggplot(results, aes(x = reorder(bee_species, p_adjusted), y = p_adjusted)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Adjusted P-values from Wilcoxon Test per Species",
    x = "Bumblebee Species",
    y = "FDR-adjusted P-value"
  ) +
  theme_minimal()


#Heatmap of Differences (F − M) per Compound
library(pheatmap)

heat_df <- paired_list %>%
  mutate(diff = F - M) %>%
  select(bee_species, compound, diff) %>%
  pivot_wider(names_from = compound, values_from = diff) %>%
  column_to_rownames("bee_species")

pheatmap(heat_df,
         main = "Sex Differences in Compound Abundance (F - M)",
         color = viridis::viridis(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE)


#Convert paired_list to long format
library(tidyverse)

# Convert to long format
paired_long <- paired_list %>%
  pivot_longer(
    cols = c(F, M),
    names_to = "bee_sex",
    values_to = "abundance"
  )

#Paired Dot Plot with Lines Connecting F ↔ M
ggplot(paired_long, aes(x = bee_sex, y = abundance, group = compound)) +
  geom_point(aes(color = bee_sex), size = 2, alpha = 0.6) +
  geom_line(aes(group = compound), alpha = 0.3) +
  facet_wrap(~ bee_species, scales = "free_y") +
  labs(
    title = "Paired Floral Scent Preferences by Sex",
    x = "Sex",
    y = "Relative Abundance"
  ) +
  theme_minimal()

#Violin + Boxplot of Abundance Distributions
ggplot(paired_long, aes(x = bee_species, y = abundance, fill = bee_sex)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  labs(
    title = "Distribution of Floral Scent Preferences by Sex",
    x = "Bumblebee Species",
    y = "Relative Abundance"
  ) +
  theme_minimal()



#option 2


# Load libraries
library(tidyverse)

# Step 1: Load the data
df <- read_csv("bee_profiles_relative_weighted_june14.csv")

# Step 2: Add a row identifier if needed
# (Assumes each row corresponds to a bee_species × sex)
df <- df %>% 
  mutate(bee_id = paste(bee_species, bee_sex, sep = "_"))

# Step 3: Pivot longer to get a long-form structure
df_long <- df %>%
  pivot_longer(
    cols = -c(bee_species, bee_sex, bee_id),
    names_to = "compound",
    values_to = "abundance"
  )


# Step: Group by species and run paired Wilcoxon test
wilcoxon_results <- df_long %>%
  pivot_wider(names_from = bee_sex, values_from = abundance) %>%
  filter(!is.na(F) & !is.na(M)) %>%
  group_by(bee_species) %>%
  summarise(
    p_value = wilcox.test(F, M, paired = TRUE)$p.value,
    statistic = wilcox.test(F, M, paired = TRUE)$statistic,
    .groups = "drop"
  )

# Step: View results
print(wilcoxon_results)



# Step 4: Pivot wider so each species/compound has separate F and M columns
paired_list <- df_long %>%
  pivot_wider(
    names_from = bee_sex,
    values_from = abundance
  ) %>%
  filter(!is.na(F) & !is.na(M))  # Keep only compound rows with both sexes present

# Step 5: View result
head(paired_list)



#option 3

# Load necessary libraries
library(tidyverse)

# Step 1: Load the dataset
df <- read_csv("bee_profiles_relative_weighted_june14.csv")

# Step 2: Create a unique ID column (optional)
df <- df %>%
  mutate(bee_id = paste(bee_species, bee_sex, sep = "_"))

# Step 3: Pivot to long format
df_long <- df %>%
  pivot_longer(
    cols = -c(bee_species, bee_sex, bee_id),
    names_to = "compound",
    values_to = "abundance"
  )

head(df_long)

write.csv(df_long, "df_long.csv")

# Step 4: Pivot to wide format to get F and M side by side
df_wide <- df_long %>%
  pivot_wider(
    names_from = bee_sex,
    values_from = abundance
  ) %>%
  filter(!is.na(F), !is.na(M))  # Only keep rows with both F and M

# Step 5: Run Wilcoxon signed-rank test for each species
wilcoxon_results <- df_wide %>%
  group_by(bee_species) %>%
  summarise(
    p_value = tryCatch(wilcox.test(F, M, paired = TRUE)$p.value, error = function(e) NA),
    statistic = tryCatch(wilcox.test(F, M, paired = TRUE)$statistic, error = function(e) NA),
    .groups = "drop"
  )

# Step 6: View results
print(wilcoxon_results)



#option 4

library(tidyverse)

# Step 1: Load the data
df <- read_csv("bee_profiles_relative_weighted_june14.csv")

# Step 2: Reshape to long format: one row per compound per bee
df_long <- df %>%
  pivot_longer(
    cols = -c(bee_species, bee_sex),
    names_to = "compound",
    values_to = "abundance"
  )

# Step 3: Create a pairing key — species + compound
# If you have a unique bee_id or sample_id, include that too
df_long <- df_long %>%
  mutate(pair_key = paste(bee_species, compound, sep = "_"))

# Step 4: Pivot to wide format: F and M side by side
df_wide <- df_long %>%
  select(pair_key, bee_sex, abundance) %>%
  pivot_wider(
    names_from = bee_sex,
    values_from = abundance
  ) %>%
  filter(!is.na(F), !is.na(M))  # keep only rows with both F and M values

# Step 5: Run the Wilcoxon signed-rank test
wilcox_result <- wilcox.test(df_wide$F, df_wide$M, paired = TRUE)

# Output the result
print(wilcox_result)


