# Supplementary Figure S5: Additional Heatmap Analysis

## Overview
This directory contains scripts for **Supplementary Figure S5**, providing complementary heatmap visualizations for chromatin state analysis across diverse samples.

## Contents

### Scripts

#### `heatmaps.R`
R script generating supplementary heatmap visualizations.
- **Purpose**: Creates detailed heatmap panels for supplement
- **Features**:
  - Multi-sample state heatmaps
  - Alternative clustering methods (hierarchical, k-means)
  - Custom color gradients
  - Annotation tracks for cell type information
  - Dendrograms showing sample relationships
- **Outputs**: Multiple heatmap variations
  - Sample-ordered heatmaps
  - Clustered sample groupings
  - State frequency heatmaps
  - Cell type × state interaction heatmaps

### Workflow

1. **Data Loading**: Import segmentation results from Figure 1
2. **Heatmap Preparation**: Format data for visualization
3. **Clustering**: Apply hierarchical clustering to identify sample groups
4. **Visualization**: Generate heatmaps with annotations
5. **Export**: Save high-resolution figures

## Analysis Details

### Heatmap Components
- **Rows**: Chromatin states (11 states)
- **Columns**: IHEC samples (154 samples)
- **Color Intensity**: Fraction or count of each state per sample
- **Annotations**: 
  - Cell type classification
  - Tissue origin
  - Quality metrics
  - Sample groupings

### Clustering Approaches
- **Hierarchical Clustering**: Distance-based grouping of similar samples
- **Distance Metrics**: Euclidean, Manhattan, Pearson correlation
- **Sample Groupings**: Identifies naturally clustering samples
- **Dendrograms**: Shows hierarchical relationships

## Usage Example

```bash
# Generate supplementary heatmaps
Rscript heatmaps.R

# Output files:
# - S5_heatmaps_clustered.pdf
# - S5_heatmaps_ordered.png
# - sample_dendrograms.pdf
```

## Requirements

- R (with ComplexHeatmap, ggplot2, tidyverse packages)
- Segmentation data (BED files or state matrices)
- Sample metadata and cell type annotations
- Color scheme definitions

## Output Files

- Clustered heatmap figures (PNG/PDF)
- Sample dendrograms (PDF)
- Cluster assignments and groupings (CSV)
- Annotation matrices (TSV/CSV)

## Biological Insights

The supplementary heatmaps reveal:
- Natural clustering of samples by cell type
- Shared chromatin state patterns within cell type groups
- Rare state combinations highlighting unique samples
- Tissue-specific state frequency patterns
- Validation of primary figure findings across diverse samples
