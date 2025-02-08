# Periodontitis Microbiome Analysis

This repository contains scripts and data for analyzing microbiome data related to periodontitis. The analysis includes data preprocessing, statistical analysis, visualization, and machine learning techniques to explore the relationship between microbial communities and periodontitis.

## Repository Structure

- **Data/**: Contains input data files (e.g., `feature-table.tsv`, `taxonomy.tsv`, `metadata_groups.tsv`).
- **Scripts/**: Contains R scripts for data analysis and visualization.
  - `boxplots_correlation_genra.R`: Generates boxplots and performs statistical tests for top genera.
  - `convert_abundances_tables.R`: Converts raw abundance data into relative abundances and merges with metadata.
  - `EDA_Periodontitis.Rmd`: Exploratory data analysis (EDA) of metadata and microbiome data.
  - `Heatmap_rel_abun_periodontitis.Rmd`: Generates heatmaps of bacterial abundances grouped by metadata categories.
  - `PCA_periodontitis.Rmd`: Performs Principal Component Analysis (PCA) on relative abundance data.
- **Results/**: Stores output files (e.g., plots, tables, and processed data).
- **comandos.txt**: Contains the pipeline commands for processing microbiome data using QIIME2.

## Data Description

- **feature-table.tsv**: ASV/OTU table with microbial abundances.
- **taxonomy.tsv**: Taxonomic classification of ASVs/OTUs.
- **metadata_groups.tsv**: Metadata file with sample information (e.g., health condition, smoking status, periodontitis extent).

## Pipeline Overview

1. **Data Preprocessing**:
   - Import raw sequencing data into QIIME2.
   - Trim primers and perform quality filtering using DADA2.
   - Generate ASV tables and representative sequences.

2. **Taxonomic Analysis**:
   - Classify ASVs using a pre-trained classifier (e.g., SILVA).
   - Filter out mitochondrial and chloroplast sequences.

3. **Statistical Analysis**:
   - Perform alpha and beta diversity analyses.
   - Generate boxplots, heatmaps, and PCA plots to explore microbial patterns.

4. **Visualization**:
   - Create barplots, heatmaps, and PCA plots to visualize results.

## Usage

### Running the Pipeline
1. Activate the QIIME2 environment:
   ```bash
   conda activate qiime2-amplicon-2024.5
