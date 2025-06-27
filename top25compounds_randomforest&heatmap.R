# Load required libraries
library(tidyverse)
library(randomForest)
library(ggplot2)

# -------------------------------
# Step 1: Load and clean the dataset
# -------------------------------
data <- read_csv("average_compound_profile_per_plant_june11.csv")

# Ensure 'Plant_species' is treated as a factor
data$Plant_species <- as.factor(data$Plant_species)

# -------------------------------
# Step 2: Prepare data
# -------------------------------

# Set plant species as the response variable
# All other columns are assumed to be compound abundances
rf_data <- data %>%
  drop_na()  # Remove any rows with NA values

# Ensure all compound columns are numeric
compound_cols <- names(rf_data)[!(names(rf_data) %in% c("Plant_species"))]
rf_data[compound_cols] <- lapply(rf_data[compound_cols], as.numeric)

# Fix column names to be valid R variable names
names(rf_data) <- make.names(names(rf_data))


# -------------------------------
# Step 3: Run Random Forest
# -------------------------------
set.seed(123)
rf_model <- randomForest(
  Plant_species ~ .,
  data = rf_data,
  importance = TRUE,
  ntree = 1000
)

# Print summary
print(rf_model)

# -------------------------------
# Step 4: Extract and plot top 20 important compounds
# -------------------------------
importance_df <- as.data.frame(importance(rf_model))
importance_df$Compound <- rownames(importance_df)

# Use MeanDecreaseGini as the importance metric
top20 <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice(1:25)

# Plot
ggplot(top20, aes(x = reorder(Compound, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 20 Indicator Compounds Across Plant Species",
    x = "Compound",
    y = "Importance (MeanDecreaseGini)"
  )


#option 2 with raw data---------------------

# Load libraries
library(tidyverse)
library(randomForest)

# Step 1: Load the dataset
data <- read_csv("plant_compound_data_raw_jun21.csv")

# Step 2: Pivot data to wide format (each row = plant individual, columns = compounds)
wide_data <- data %>%
  pivot_wider(
    id_cols = c(Individual, Plant_species),
    names_from = Compound,
    values_from = Abundance,
    values_fill = 0
  )

wide_data <- data %>%
  pivot_wider(
    id_cols = c(Individual, Plant_species),
    names_from = Compound,
    values_from = Abundance,
    values_fill = list(Abundance = 0)
  )

pivot_wider(names_from = Compound, values_from = Abundance, values_fill = list(Abundance = 0))


head(data)

# Step 3: Prepare data for Random Forest
rf_data <- wide_data %>% select(-Individual)
rf_data$Plant_species <- as.factor(rf_data$Plant_species)

# Step 4: Run Random Forest classification
set.seed(123)
rf_model <- randomForest(
  Plant_species ~ .,
  data = rf_data,
  importance = TRUE,
  ntree = 1000
)

# Step 5: View model output
print(rf_model)

# Step 6: Extract variable importance
importance_df <- as.data.frame(importance(rf_model))
importance_df$Compound <- rownames(importance_df)

# Step 7: Plot top 20 compounds by MeanDecreaseGini
top_compounds <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice(1:20)

ggplot(top_compounds, aes(x = reorder(Compound, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 20 Indicator Compounds Across Plant Species",
    x = "Compound",
    y = "Importance (MeanDecreaseGini)"
  )

#Heatmap 25 most important compounds

# Load necessary libraries
library(tidyverse)
library(pheatmap)
library(viridis)

# Step 1: Load plant compound data
compound_data <- read_csv("average_compound_profile_per_plant_june11.csv")

# Step 2: Set plant species as rownames
compound_data <- compound_data %>%
  mutate(Plant_species = gsub(" ", "_", Plant_species)) %>%
  column_to_rownames("Plant_species")

# Step 3: Load or reuse top 20 compounds (from previous Random Forest analysis)
# If you already saved importance_df from the previous RF model:
# top_compounds <- importance_df %>%
#   arrange(desc(MeanDecreaseGini)) %>%
#   slice(1:20) %>%
#   pull(Compound)

# For now, assume you manually define the vector:
top_compounds <- c("D-Limonene_B_1028", "Phenylethyl Alcohol_B_1110", 
                   "Benzaldehyde, 2,5-dimethyl-_1235", "Naphthalene_1191",
                   "1-Octen-3-ol_978", "Benzaldehyde, 2,4-dimethyl-_B_1218",
                   "Mesitylene_A_992", "Butanoic acid, butyl ester_994",
                   "Octanal_1002", "Phenyl-β-D-glucoside_974", "Acetophenone_1063",
                   "Benzaldehyde_B_959", "Benzene, 1-methyl-2-propyl-_1049",
                   "1-Nonanol, 4,8-dimethyl-_1069", "α-Pinene_932", "Geranyl acetone_B_1445",
                   "o-Cymene_1023", "3-Carene_A_1006", "Farnesan_C_1375",
                   "Farnesan_A_1051", "cis-β-Farnesene_A_1452", 
                   "4-Hexen-1-ol, acetate-B_1003", "Linalyl acetate_A_1034",
                   "Geranyl isovalerate_B_1462", "Mesitylene_B_1018")

# Step 4: Subset the data to top 20 compounds
compound_matrix <- compound_data[, top_compounds]

compound_matrix_mat <- as.matrix(compound_matrix)

write.csv(compound_matrix, "25important_compounds.csv")

# Step 5: Create the heatmap
pheatmap(
  mat = compound_matrix,
  color = plasma(99),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",  # standardize rows for better contrast
  fontsize_row = 9,
  fontsize_col = 9,
  border_color = NA,
  main = "Heatmap of Top 20 VOCs Across Plant Species"
)

#heatmap adjusted

# Load required libraries
library(tidyverse)
library(pheatmap)
library(viridis)

# Step 1: Load data
compound_matrix <- read_csv("top25_compounds.csv")  # Replace with your file if different
compound_matrix <- column_to_rownames(compound_matrix, var = "Plant_species")

# Step 2: Reorder plant species
plant_species_order <- c(
  "Ranunculus acris", "Tilia sp2", "Epilobium angustifolium", "Veronica spicata",
  "Rhinanthus minor", "Nepeta racemosa", "Lavandula angustifolia", "Stachys palustris",
  "Ballota nigra", "Buddleja davidii", "Ligustrum vulgare", "Symphoricarpos albus",
  "Hydrangea paniculata", "Heracleum sphondylium", "Jacobaea vulgaris", "Centaurea nigra",
  "Cirsium vulgare", "Cirsium arvense", "Arctium minus", "Rubus fruticosus agg.",
  "Rubus caesius", "Trifolium repens", "Lotus corniculatus", "Hypericum perforatum",
  "Geranium pratense", "Impatiens glandulifera", "Cornus sanguinea", "Calystegia silvatica"
)

compound_matrix <- compound_matrix[plant_species_order, ]

# Step 3: Create custom color palette with zero color
zero_color <- "#fffafe"  # Custom color for zero values
color_palette <- plasma(99)
custom_colors_final <- c(zero_color, color_palette)

# Step 4: Define color breaks (include zero break)
breaks <- c(0, seq(0.0001, 0.5, length.out = 100))

# Step 5: Plot heatmap
pheatmap(
  mat = compound_matrix,
  color = custom_colors_final,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 90,  # Rotate compound names
  border_color = "grey90",
  main = "Relative Abundance of Top 20 Indicator Compounds"
)
