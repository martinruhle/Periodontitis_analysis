library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# Load feature table (OTU/ASV table)
feature_table <- read.delim("C:/Users/User/Documents/Ecología_microbiana/Proyecto_final/GitHub_repository/Data/feature-table.tsv", row.names = 1)

# Load taxonomy data
taxonomy <- read.delim("C:/Users/User/Documents/Ecología_microbiana/Proyecto_final/GitHub_repository/Data/taxonomy.tsv", header = TRUE, sep = "\t")

# Load metadata file
metadata <- read.csv("C:/Users/User/Documents/Ecología_microbiana/Proyecto_final/GitHub_repository/Data/metadata_groups.tsv", sep = "\t")

# FIltering metdata variables
metadata_fil <- metadata %>%
  select(sample.id, AGE, BOP_total, CAL_total, Cigarrettes_day, health_condition, Months_smoking, geo_loc_name_country, health_condition, RUN_ID, sex, Smoking, teeth, Periodontitis_extent_calc, Periodontitis_extent, Periodontitis_severity, cigarrettes_group, Months_smoking_group, teeth_group, BOP_group, CAL_group)

# Merge feature table with taxonomy
merged_table <- feature_table %>%
  rownames_to_column("Feature.ID") %>%
  left_join(taxonomy, by = c("Feature.ID" = "Feature.ID"))

# Specify the columns to move (126 to 129)
cols_to_move <- names(merged_table)[126:127]

# Specify the new order: 1 (keep the first column), then the columns to move, then the rest
new_order <- c(names(merged_table)[1], cols_to_move, names(merged_table)[-c(1, 126:127)])

# Reorder the data frame
merged_table <- merged_table[, new_order]

# Split the "Taxon" column into separate columns
merged_table <- merged_table %>%
  separate(
    col = Taxon,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specie"),
    sep = ";",
    fill = "right",   # Ensures missing values are filled with NA
    extra = "warn"    # Warns if there are extra levels in the Taxon
  )

# Optional: Remove prefixes like "d__", "p__", etc.
merged_table <- merged_table %>%
  mutate(across(
    c(Domain, Phylum, Class, Order, Family, Genus, Specie),
    ~ sub("^[a-z]__*", "", .)
  ))

merged_table_conf <- merged_table %>%
  filter(Confidence > 0.95 & !is.na(Genus))

# Transform the table to long format, keeping only relevant columns
long_table <- merged_table_conf %>%
  select(Genus, starts_with("SRR")) %>%
  pivot_longer(
    cols = starts_with("SRR"),
    names_to = "sample.id",
    values_to = "Abundance"
  )

# Summarize to ensure unique Sample-Genus pairs
long_table <- long_table %>%
  group_by(sample.id, Genus) %>%
  summarize(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

# Transform the table to wide format
wide_table <- long_table %>%
  pivot_wider(
    names_from = Genus,
    values_from = Abundance,
    values_fill = list(Abundance = 0) # Specify fill for Abundance explicitly
  )

# Convert to relative abundances
relative_table <- wide_table %>%
  mutate(across(-sample.id, ~ . / rowSums(across(-sample.id))))

# Join the datasets
rel_abun_metadata <- metadata_fil %>%
  inner_join(relative_table, by = "sample.id")

write.csv(x = rel_abun_metadata, file = "C:/Users/User/Documents/Ecología_microbiana/Proyecto_final/GitHub_repository/Data/genus_rel_abund_meta.csv")

# Transform the dataset to long format, specifying columns 8 onward for species
rel_abun_long <- rel_abun_metadata %>%
  pivot_longer(cols = 21:ncol(rel_abun_metadata), names_to = "Genus", values_to = "abundance")


