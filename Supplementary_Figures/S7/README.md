# Supplementary Figure S7: HMM Emission Matrix Analysis

## Overview
This directory contains scripts for **Supplementary Figure S7**, focusing on hidden Markov model parameters and emission probability analysis.

## Contents

### Scripts

#### `emission_matrix.r`
R script analyzing HMM emission matrices and state-signal relationships.
- **Purpose**: Visualizes how chromatin states relate to observed histone signals
- **Analysis**:
  - Extracts emission probability matrices from trained HMM models
  - Creates publication-quality heatmaps
  - Quantifies state-signal specificity
  - Identifies signature marks for each state
  - Compares emission profiles across models
- **Visualizations**: Heatmaps showing state × histone mark relationships
- **Output**: Emission probability figures, state signature plots

## Analysis Details

### Emission Probability Concept
- **Definition**: P(histone mark signal | chromatin state)
- **Interpretation**: 
  - High probability = mark is typically present in that state
  - Low probability = mark is typically absent in that state
  - Specific patterns define state biological identity
- **Application**: Validates that states have distinguishing signal profiles

### Heatmap Components
- **Rows**: Chromatin states (11 states)
- **Columns**: Histone marks (6 marks: H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3, H3K9me3)
- **Color Intensity**: Emission probability (0 to 1 scale)
- **Annotations**: State names, biological classifications

### State Signatures
- **Active Promoters**: High H3K4me3, moderate H3K27ac
- **Active Enhancers**: High H3K27ac, moderate H3K4me1
- **Repressive**: High H3K27me3 or H3K9me3
- **Poised States**: H3K4me3 + H3K27me3 bivalence
- **Heterochromatin**: High H3K9me3
- **Background**: Low signal across all marks

## Usage Example

```bash
# Generate emission matrix analysis
Rscript emission_matrix.r

# Output files:
# - S7_emission_matrix_chromhmm.pdf
# - S7_emission_matrix_episegmix.pdf
# - S7_emission_matrix_episegmix_meth.pdf
# - state_signatures.csv
```

## Requirements

- R (with ComplexHeatmap, ggplot2, tidyverse packages)
- HMM model files with emission parameters
- Model training data or pre-computed emission matrices
- State definitions and naming schemes

## Data Inputs

- Trained HMM parameters from each model:
  - ChromHMM 10-state model
  - EpiSegMix 10-state model
  - EpiSegMix with methylation model

## Output Files

- Emission probability heatmaps (PNG/PDF)
  - One per model (ChromHMM, EpiSegMix, EpiSegMix+Meth)
  - High resolution for publication
- State signature summaries (TSV)
  - Top histone marks per state
  - Emission probability scores
- Model comparison plots
  - Overlaid emission profiles
  - Divergence identification

## Biological Validation

The emission matrix analysis confirms:
- States have biologically meaningful chromatin signatures
- Distinct histone mark combinations define state identity
- Model learned patterns consistent with known chromatin biology
- Methylation integration captures additional information
- State classification is reproducible and robust

## Interpretation Guide

**High Emission Values Indicate:**
- Frequent co-occurrence of state and signal
- State-specific chromatin feature
- Reliable state identification marker

**Low Emission Values Indicate:**
- Rare co-occurrence
- Non-defining feature for that state
- May be present in minority of state instances

**Across-model Comparisons:**
- Similar emission patterns validate consistency
- Divergent patterns highlight model-specific features
- Methylation data impact on EpiSegMix+Meth variant
