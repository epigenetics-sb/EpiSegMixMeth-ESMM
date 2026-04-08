# Figure 2: State Overlap and Heatmap Analysis

## Overview
This directory contains scripts for **Figure 2**, which visualizes the overlap between different segmentation models (ChromHMM and EpiSegMix) and generates heatmaps showing state concordance across IHEC samples.

## Contents

### Scripts

#### `state_overlap.ns.r`
R script analyzing state overlap between ChromHMM and EpiSegMix models.
- **Purpose**: Quantifies concordance and divergence between segmentation outputs
- **Analysis**: Calculates overlap statistics using contingency tables
- **Output**: State-by-state comparison matrices
- **Features**: No-signal state handling, normalization options

#### `ChromHMM_OverlapNS.R`
R script generating overlap visualizations with no-signal (NS) state filtering.
- **Purpose**: Creates publication-quality overlap plots
- **Filtering**: Excludes background/no-signal states for cleaner visualization
- **Visualization**: Heatmaps, bar plots, correlation matrices
- **Output**: PNG and PDF formatted figures

#### `heatmaps.R`
Comprehensive heatmap generation script for state comparisons.
- **Purpose**: Creates multi-sample heatmaps showing state patterns
- **Features**: 
  - Color-coded state representations
  - Sample grouping and sorting
  - Hierarchical clustering
  - Scale normalization
- **Output**: High-resolution heatmap figures

#### `pygenome_commands.sh`
Shell script with commands for pyGenome-based heatmap generation.
- **Purpose**: Alternative visualization approach using pyGenome tools
- **Usage**: Command reference for genomic annotation visualization

### Data Files

#### Overlap Analysis Results
- **EpiSegMixMeth_overlap_ChromHMM-NS.states.txt**: State overlap data with no-signal filtering
- **EpiSegMixMeth_ChromHMM_NS-overlap_counts.csv**: Summary statistics of state overlaps

#### Supporting Files
- **States_tables_with colors.xlsx**: Excel reference file with state definitions and color schemes
  - Contains standardized color palettes for consistent visualization
  - State nomenclature and biological meanings
  - Cell type-specific state enrichments

### Workflow

1. **Data Input**: Load segmentation outputs from Figure 1 analysis
2. **Overlap Calculation**: Run `state_overlap.ns.r` to compute concordance metrics
3. **Visualization**: Execute `ChromHMM_OverlapNS.R` to generate publication figures
4. **Heatmap Generation**: Use `heatmaps.R` for detailed state pattern visualization
5. **Quality Control**: Verify output figures and statistics

## Analysis Details

### State Overlap Metrics
- **Concordance Score**: Fraction of genome with identical state assignments
- **Jaccard Index**: Measures similarity between state assignments
- **Confusion Matrix**: Detailed comparison of state assignments between models
- **No-Signal Handling**: Special treatment for background states to avoid bias

### Heatmap Features
- **Sample-by-State Matrix**: Rows = samples, columns = chromatin states
- **Color Intensity**: Represents abundance or enrichment of each state
- **Hierarchical Clustering**: Identifies similar samples and state patterns
- **Cell Type Annotations**: Shows grouping by cell type or tissue

## Usage Example

```bash
# Calculate overlap statistics
Rscript state_overlap.ns.r

# Generate overlap visualization
Rscript ChromHMM_OverlapNS.R

# Create comprehensive heatmaps
Rscript heatmaps.R

# Optional: Generate pyGenome visualizations
bash pygenome_commands.sh
```

## Requirements

- R (with ggplot2, reshape2, ComplexHeatmap, tidyverse packages)
- Input CSV/TXT files from Figure 1 classification
- Color scheme definition file (States_tables_with colors.xlsx)
- PyGenome (optional, for alternative visualizations)

## Output Files

- State overlap matrices (CSV format)
- Heatmap visualizations (PNG/PDF)
- Concordance statistics plots
- Color-coded state comparison figures
- Manuscript-ready figure panels

## Key Findings

The analysis demonstrates:
- High overall concordance between EpiSegMix and ChromHMM models
- Specific regions of model divergence related to weak chromatin signals
- State-specific enrichments in particular cell types
- Importance of methylation data for accurate segmentation in certain contexts
