# Supplementary Figure S10: PyGenome Extended Visualizations

## Overview
Directory for **Supplementary Figure S10**, containing pyGenome visualization scripts and configuration files.

## Contents

### Scripts

#### `pygenome_commands.sh`
Extended pyGenome commands for alternative visualization approaches.
- **Purpose**: Generates epilogos plots and custom genomic visualizations
- **Features**: Custom annotation tracks, multi-sample views
- **Output**: High-resolution genomic landscape plots

## PyGenome Integration

PyGenome provides:
- **Epilogos Plots**: Multi-sample chromatin state visualization
- **Landscape Views**: Horizontal genome browser-like representation
- **Feature Overlays**: Integration with genomic annotations
- **Publication Quality**: High-resolution output suitable for journals

## Analysis Features

- Whole-genome visualization across all samples
- State frequency distributions by genomic position
- Integrated gene annotation overlays
- Enrichment highlighting of specific regions
- Color-coded by state type for clarity

## Usage

```bash
# Run pyGenome visualizations
bash pygenome_commands.sh

# Generated output:
# - Epilogos plots (PNG)
# - Landscape views (PDF/PNG)
# - Combined visualization panels
```

## Requirements

- PyGenome installation and dependencies
- Segmentation BED files
- Gene annotation files
- Reference genome sequence
- Color configuration files

## Output Files

- Epilogos plots across all samples (PNG)
- Landscape-view segmentations (PDF)
- Feature overlay plots
- State distribution visualizations
- Manuscript-ready figure components

## Key Advantages

- **Comprehensive View**: All samples and regions in single visualization
- **Publication Quality**: High resolution, professional formatting
- **Interpretability**: Clear state color coding and annotations
- **Biological Context**: Integration with genomic features
- **Reproducibility**: Scripted generation ensures consistency
