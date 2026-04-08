# Supplementary Figures S6-S8: Advanced Statistical Analysis

## Overview
This directory contains scripts for **Supplementary Figures S6-S8**, which present additional statistical analyses including Jaccard indices, emission matrices, and specialized filtering approaches.

## S6 Contents

### Scripts

#### `jaccard.r`
R script focusing on Jaccard index calculations and visualizations.
- **Purpose**: Detailed analysis of state similarity between models
- **Calculations**:
  - Pairwise Jaccard indices
  - Matrix visualization
  - Distribution analysis
  - State-specific comparisons
- **Output**: Jaccard heatmaps, distribution plots

#### `JS_distance_and_coverage.r`
R script analyzing Jensen-Shannon distances and genomic coverage.
- **Purpose**: Alternative similarity metric and coverage analysis
- **Analysis**:
  - JS divergence between segmentation models
  - Genomic coverage statistics
  - State-specific coverage profiles
  - Comparative coverage visualization
- **Metrics**:
  - JS Distance: Information-theoretic similarity measure
  - Coverage: Base pair counts per state
  - Enrichment: Relative abundance across models
- **Output**: Coverage comparison plots, JS distance matrices

### S7 Contents

#### `emission_matrix.r`
R script analyzing HMM emission matrices and state probabilities.
- **Purpose**: Examines hidden Markov model parameters
- **Analysis**:
  - Emission probability distributions
  - State-to-observation mappings
  - Transition probabilities
  - Model parameter validation
- **Visualization**: Heatmaps of emission probabilities
- **Output**: Emission matrix figures, parameter summary tables

### S8 Contents

#### `circlize_plot.r`
R script creating circular/circos plots for chromatin state data.
- **Purpose**: Alternative visualization format for state distributions
- **Features**:
  - Circular genome representation
  - State coloring by chromosomal position
  - Multiple track display
  - Chord diagrams for state transitions
  - Interactive HTML output
- **Output**: Circos plots (PNG), interactive circos (HTML)

## Workflow

1. **S6 - Jaccard Analysis**: Run `jaccard.r` for model comparison
2. **S6 - Coverage Analysis**: Execute `JS_distance_and_coverage.r`
3. **S7 - Emission Analysis**: Run `emission_matrix.r` for model parameters
4. **S8 - Circular Visualization**: Execute `circlize_plot.r` for genome plots

## Analysis Details

### S6: Jaccard and Coverage Analysis
- **Jaccard Index Range**: 0 (no overlap) to 1 (perfect agreement)
- **State-Specific Jaccard**: Computed independently for each chromatin state
- **Coverage Metrics**:
  - Genomic coverage (percentage of genome per state)
  - State frequencies across samples
  - Cell type-specific coverage patterns

### S7: Emission Matrix Analysis
- **HMM Framework**: Models chromatin state transitions
- **Emission Probabilities**: 
  - Probability of observing histone mark signals given a state
  - Defines state-specific histone mark profiles
  - Parameter validation against observed data
- **Key Parameters**: Transition probabilities, emission distributions

### S8: Circular Visualization
- **Circos Components**:
  - Outer ring: Genomic coordinates and chromosomes
  - Inner tracks: Chromatin states by chromosome
  - Chord connections: State transitions between regions
  - Color coding: By state category or model
- **Advantages**: 
  - Whole-genome visualization
  - Pattern recognition across all chromosomes
  - Intuitive representation of state arrangements

## Usage Example

```bash
# S6: Jaccard analysis
Rscript jaccard.r

# S6: Coverage and JS distance
Rscript JS_distance_and_coverage.r

# S7: Emission matrix analysis
Rscript emission_matrix.r

# S8: Circular plots
Rscript circlize_plot.r
```

## Requirements

- R (with tidyverse, ggplot2, circlize, igraph packages)
- Segmentation BED files
- HMM model parameters and emission matrices
- Reference genome coordinates (chromosome sizes)
- State definitions and color schemes

## Data Inputs

- State overlap matrices (from Figure 2 analysis)
- HMM parameters and emission probability matrices
- Genomic coverage files (BigWig format optional)
- Segmentation outputs from all three models

## Output Files

### S6
- Jaccard index heatmaps (PNG/PDF)
- Jaccard distribution plots
- JS distance matrices (CSV/TSV)
- Coverage comparison figures

### S7
- Emission probability heatmaps (PNG/PDF)
- Parameter validation plots
- Model specification summaries

### S8
- Circos plots - full genome (PNG)
- Interactive circos visualizations (HTML)
- State distribution by chromosome plots
- Chord diagrams showing state transitions

## Key Insights

The supplementary analyses demonstrate:
- Strong Jaccard indices confirm model agreement
- JS divergence identifies specific regions of model divergence
- HMM emission matrices validate state-specific histone mark profiles
- Circular plots reveal chromosomal distribution patterns
- Consistent state signatures across genomic regions
