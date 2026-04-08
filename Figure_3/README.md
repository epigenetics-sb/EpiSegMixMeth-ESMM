# Figure 3: Genomic Coverage, Jaccard Index, and Gene Expression Analysis

## Overview
This directory contains scripts for **Figure 3**, which analyzes genomic coverage patterns, computes Jaccard similarity indices between models, and correlates segmentation states with gene expression data.

## Contents

### Scripts

#### `genomics_coverage_with_jaccard.r`
Comprehensive R script analyzing genomic coverage and Jaccard index calculations.
- **Purpose**: Quantifies similarity between segmentation models using Jaccard index
- **Analysis**:
  - Calculates Jaccard similarity for each chromatin state
  - Compares EpiSegMix vs ChromHMM genome-wide coverage
  - Generates coverage distribution plots
  - Creates Jaccard heatmaps and bar plots
- **Visualization**: Publication-quality figures showing coverage statistics
- **Output**: Jaccard matrices, coverage comparison plots

#### `gene_expression_total-RNA.r`
R script correlating segmentation states with gene expression data.
- **Purpose**: Links chromatin states to gene activity levels
- **Analysis**:
  - Assigns genes to chromatin states based on overlap
  - Computes average expression per state
  - Correlates state predictions with RNA-seq measurements
  - Generates predictive accuracy plots
- **Data Sources**: 
  - Total RNA-seq data (200bp resolution)
  - State predictions from all three models (ChromHMM, EpiSegMix, EpiSegMix+Meth)
- **Output**: Gene expression vs chromatin state plots, prediction accuracy metrics

#### `heatmap_sankey.r`
R script generating Sankey diagrams and heatmaps for state transitions.
- **Purpose**: Visualizes relationships and transitions between chromatin states
- **Features**:
  - Sankey flow diagrams showing state transitions
  - Hierarchical heatmaps of state relationships
  - Interactive visualizations (HTML widget format)
  - Color-coded flow by state type
- **Output**: Interactive HTML visualizations, static PNG/PDF figures

#### `jaccard_index_caclulation.sh`
Bash script automating Jaccard index calculations across genomic regions.
- **Purpose**: Parallel computation of Jaccard indices
- **Features**: Batch processing, memory-efficient region analysis
- **Output**: Jaccard index results for each region/state

### Workflow

1. **Genomic Coverage Analysis**: Run `genomics_coverage_with_jaccard.r` for initial comparison
2. **Gene Expression Correlation**: Execute `gene_expression_total-RNA.r` to link states to expression
3. **State Relationship Visualization**: Use `heatmap_sankey.r` for flow visualization
4. **Index Calculations**: Automated computation via `jaccard_index_caclulation.sh`

## Analysis Details

### Jaccard Index Calculation
- **Definition**: J(A,B) = |A ∩ B| / |A ∪ B|
- **Interpretation**: 
  - Value 1.0 = perfect agreement between models
  - Value 0.0 = no overlap between models
- **State-Specific Analysis**: Computed separately for each chromatin state
- **Genome-Wide Summary**: Overall model concordance metrics

### Gene Expression Analysis
- **Data Resolution**: 200bp bins for expression prediction
- **Models Compared**:
  - ChromHMM 10-state predictions
  - EpiSegMix 10-state predictions  
  - EpiSegMix with methylation predictions
- **Prediction Accuracy**: Evaluated using correlation metrics
- **Cell Type Specificity**: Expression patterns across different cell types

### Sankey Diagram Features
- **Node Size**: Represents frequency or importance of each state
- **Flow Width**: Indicates strength of state transitions or relationships
- **Color Coding**: Distinct colors for different state categories
- **Interactive Format**: HTML widgets allow zoom, pan, and hover details

## Usage Example

```bash
# Analyze coverage and Jaccard indices
Rscript genomics_coverage_with_jaccard.r

# Correlate segmentation with gene expression
Rscript gene_expression_total-RNA.r

# Generate state transition visualizations
Rscript heatmap_sankey.r

# Optional: Parallel Jaccard calculations
bash jaccard_index_caclulation.sh
```

## Requirements

- R (with ggplot2, reshape2, ggalluvial, tidyverse packages)
- Input BED files with segmentation states
- RNA-seq data in standard format (CPM, RPKM, or raw counts)
- Gene annotation file (GTF/GFF format recommended)
- Bash shell for parallel processing

## Data Inputs

- Segmentation files (BED format) from Figure 1
- RNA-seq expression matrix (samples × genes)
- Gene coordinate annotations
- Chromatin state definitions and color schemes

## Output Files

- Jaccard index matrices (CSV)
- Coverage comparison plots (PNG/PDF)
- Gene expression correlation plots
- Sankey diagrams (HTML interactive, PNG static)
- Prediction accuracy metrics
- State-specific gene expression heatmaps
- Manuscript-ready figure panels

## Key Findings

The analysis reveals:
- High Jaccard index scores between models, indicating substantial agreement
- Significant correlation between chromatin states and gene expression levels
- Enhanced prediction accuracy with methylation-aware segmentation
- Cell type-specific state-expression relationships
- Distinct chromatin state profiles associated with different gene activity levels
