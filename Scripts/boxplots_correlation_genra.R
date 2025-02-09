library(ggpubr)

# rel_abun_long generated in "convert_abundances_table.R"

# Identify top 10 genera by mean abundance
top_10_genus <- rel_abun_long %>%
  group_by(Genus) %>%
  summarise(total_mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(order_by = total_mean_abundance, n = 10) %>%
  pull(Genus)

# List of categorical variables
categorical_vars <- c("Periodontitis_extent", "cigarrettes_group", "Months_smoking_group", "BOP_group", "CAL_group")

######

# Loop through categorical variables
for (var in categorical_vars) {
  
  # Filter for top 10 genera
  data_subset <- rel_abun_long %>%
    filter(Genus %in% top_10_genus)
  
  # Calculate pairwise Wilcoxon tests and prepare labels for plotting
  test_results <- data_subset %>%
    group_by(Genus) %>%
    do({
      pairwise_results <- pairwise.wilcox.test(
        .$abundance, .[[var]], # Fixed syntax
        p.adjust.method = "bonferroni"
      )
      expand.grid(
        group1 = rownames(pairwise_results$p.value),
        group2 = colnames(pairwise_results$p.value)
      ) %>%
        mutate(
          p_value = as.vector(pairwise_results$p.value),
          significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01  ~ "**",
            p_value < 0.05  ~ "*",
            TRUE            ~ "ns"
          ),
          Genus = .$Genus[1] # Add genus information
        )
    }) %>%
    filter(!is.na(p_value))
  
  # Add significance labels to the boxplots
  for (genus in unique(test_results$Genus)) {
    genus_data <- data_subset %>% filter(Genus == genus)
    genus_test <- test_results %>% filter(Genus == genus)
    
    # Boxplot with significance annotations
    boxplot <- genus_data %>%
      ggplot(aes(x = !!sym(var), y = abundance, fill = !!sym(var))) +
      geom_boxplot() +
      geom_text(
        data = genus_test,
        aes(
          x = match(group1, levels(genus_data[[var]])), # Map group1 to x
          y = max(genus_data$abundance) * 1.1,          # Adjust for label placement
          label = significance
        ),
        inherit.aes = FALSE,
        vjust = -0.5
      ) +
      labs(
        title = paste("Boxplot for", genus, "with Pairwise Significance"),
        x = var, y = "Abundance"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(boxplot)
  }
}

######

# Create a list to store matrices for each genus
genus_matrices <- list()

# Loop through each genus
for (genus in top_10_genus) {
  
  # Filter data for the current genus
  genus_data <- rel_abun_long %>% filter(Genus == genus)
  
  # Initialize a list to store matrices for each categorical variable
  genus_matrices[[genus]] <- list()
  
  # Loop through each categorical variable
  for (var in categorical_vars) {
    
    # Perform pairwise Wilcoxon test
    pairwise_results <- pairwise.wilcox.test(
      genus_data$abundance, genus_data[[var]], 
      p.adjust.method = "bonferroni"
    )
    
    # Create a matrix from the p-values
    p_value_matrix <- pairwise_results$p.value
    p_value_matrix[is.na(p_value_matrix)] <- 1  # Replace NA with 1 for clarity
    
    # Store the matrix in the list
    genus_matrices[[genus]][[var]] <- p_value_matrix
  }
}

# Create a list to store matrices for each species
species_matrices <- list()

# Loop through the top 10 species
for (species in top_10_genus) {
  
  # Filter data for the current species
  species_data <- rel_abun_long %>% filter(Genus == species)
  
  # Initialize a list to store matrices for each categorical variable
  species_matrices[[species]] <- list()
  
  # Loop through all categorical variables
  for (var in categorical_vars) {
    
    # Perform pairwise Wilcoxon test
    pairwise_results <- pairwise.wilcox.test(
      species_data$abundance, species_data[[var]], 
      p.adjust.method = "bonferroni"
    )
    
    # Create a matrix from the p-values
    p_value_matrix <- pairwise_results$p.value
    p_value_matrix[is.na(p_value_matrix)] <- 1  # Replace NA with 1 for clarity
    
    # Store the matrix in the list
    species_matrices[[species]][[var]] <- p_value_matrix
  }
}

# Loop through the species and their respective categorical variable matrices
for (species in names(species_matrices)) {
  cat("\n==========", species, "==========\n") # Header for the species
  
  for (var in names(species_matrices[[species]])) {
    cat("\n---", var, "---\n") # Sub-header for the categorical variable
    
    # Print the p-value matrix
    print(species_matrices[[species]][[var]])
  }
}
