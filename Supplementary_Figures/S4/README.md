# Supplementary Figures S4-S5: State Distribution and Heatmap Variations

## Overview
This directory contains scripts and data for **Supplementary Figures S4-S5**, which provide additional heatmap visualizations and state distribution analyses complementing the main figures.

## S4 Contents

### Scripts

#### `merge_rna_ini.sh`
Bash script for merging RNA-seq data with INI configuration files.
- **Purpose**: Integrates gene expression data with segmentation configurations
- **Operations**: Data concatenation, format standardization, validation

#### `pygenome_commands.sh`
Shell script with pyGenome commands for figure generation.
- **Purpose**: Command reference for genomic visualization tools

### S5 Contents

#### `heatmaps.R`
R script generating supplementary heatmap variations.
- **Purpose**: Creates alternative heatmap visualizations
- **Features**:
  - Different clustering approaches
  - Alternative color schemes
  - Subset-specific visualizations
  - High-resolution outputs
- **Output**: High-quality heatmap figures

## Workflow

1. **Data Integration**: Merge RNA and segmentation data via shell script
2. **Heatmap Generation**: Generate variations using `heatmaps.R`
3. **Figure Export**: Create publication-quality PNG and PDF outputs

## Requirements

- Bash shell
- R (with ComplexHeatmap, tidyverse packages)
- RNA-seq data files
- INI configuration files from Figure 1

## Output Files

- Alternative heatmap visualizations (PNG/PDF)
- State distribution supplementary figures
- Integrated data matrices (CSV)
- Clustered sample groupings

## Key Content

- Additional heatmap perspectives on state distributions
- Alternative clustering and ordering schemes
- Supplementary validation of main figure patterns
- Extended cell type and state relationships
