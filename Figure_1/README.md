# Figure 1: EpiSegMix Classification and Random Forest Analysis

## Overview
This directory contains scripts and data for **Figure 1** of the manuscript, which demonstrates the classification accuracy and feature importance of the EpiSegMix segmentation model compared to ChromHMM using Random Forest machine learning approaches.

## Contents

### Scripts

#### `classifier.sh`
Main bash script for running the EpiSegMix and ChromHMM classifiers on input segmentation files.
- **Purpose**: Automates the classification pipeline for genomic segmentation data
- **Parameters**: Accepts BED files, state numbers, output directories, histone marks, and sample names
- **Dependencies**: R scripts, reference files, singularity containers
- **Output**: Classification results in CSV format

#### `classifier_slurm.sh`
SLURM job submission wrapper for the classifier.sh script with array job support.
- **Purpose**: Enables parallel processing of 154 IHEC samples on HPC clusters
- **Features**: Array job parallelization, memory management, error handling
- **Job Configuration**: 16 CPUs, 40GB RAM per task

#### `relable.r`
R script for relabeling and processing classifier results.
- **Purpose**: Post-processes classification outputs and handles state renaming
- **Functions**: Data aggregation, state consolidation

#### `relable_slurm.sh`
SLURM wrapper for relabeling operations with array job support.
- **Purpose**: Parallel execution of relabeling operations across samples
- **Usage**: Submits array jobs for efficient processing

#### `Random_forest_classifier.R`
Comprehensive R script implementing Random Forest classification analysis.
- **Purpose**: Trains and evaluates Random Forest models for EpiSegMix vs ChromHMM classification
- **Models**: Implements three classification models:
  - **CHMM**: Standard ChromHMM 10-state model
  - **ESMM**: EpiSegMix 10-state model
  - **ESMM_METH**: EpiSegMix with methylation data (WGBS)
- **Analysis**: Feature importance visualization, model comparison, accuracy metrics
- **Output**: Classification accuracy plots, feature importance rankings

#### `154_samples.R`
R script analyzing classification results across all 154 IHEC samples.
- **Purpose**: Aggregates and visualizes classification performance statistics
- **Output**: Sample-level classification accuracy summaries

#### `pygenome_commands.sh`
Shell script with pyGenome package commands for visualization.
- **Purpose**: Generates epilogos plots and genomic visualizations
- **Usage**: Command reference for pyGenome-based figure generation

### Data Files

#### Configuration Files
- **config.txt**: Configuration parameters for batch processing
- **episegmixmeth.ini**: EpiSegMix with methylation pipeline configuration
- **episegmixmeth_relabelled.ini**: Post-processing configuration
- **H3K*.signals.bigwig.ini**: Histone mark signal files (H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3, H3K9me3)
- **WGBS.signals.bigwig.ini**: Whole-genome bisulfite sequencing methylation signals

#### Classification Results
- **IHEC_ChromHMM_10N_all_states_samples_classifier.csv**: ChromHMM classification results (154 samples)
- **IHEC_EpiSegMix_10N_all_states_samples_classifier.csv**: EpiSegMix classification results
- **IHEC_EpiSegMix_meth_10N_all_states_samples_classifier.csv**: EpiSegMix+methylation results

#### Annotation Files
- **DEEP_Bcells_trained_*.csv**: DEEP cell line training annotations (ChromHMM, EpiSegMix, EpiSegMix_meth)
- **DEEP_samples_annotation_summary.tab**: DEEP samples summary table

### Workflow

1. **Data Preparation**: Prepare 154 IHEC sample segmentation files in BED format
2. **Classification**: Run `classifier_slurm.sh` to classify all samples (parallel processing)
3. **Relabeling**: Execute `relable_slurm.sh` to standardize state labels
4. **Analysis**: Run `Random_forest_classifier.R` to train models and generate figure panels
5. **Visualization**: Use `154_samples.R` and `pygenome_commands.sh` for comprehensive visualization

## Key Features

- **Multi-model Comparison**: Compares three distinct segmentation approaches
- **Feature Importance**: Analyzes which histone marks contribute most to classification
- **Parallel Processing**: SLURM integration for processing 154 samples efficiently
- **Methylation Integration**: Includes methylation-aware EpiSegMix variant
- **Comprehensive Metrics**: Accuracy, precision, recall, and feature rankings

## Usage Example

```bash
# Submit classification jobs for all 154 samples
sbatch classifier_slurm.sh

# After classification completes, submit relabeling jobs
sbatch relable_slurm.sh

# Generate figure panels
Rscript Random_forest_classifier.R
```

## Requirements

- R (with randomForest, ggplot2, reshape2, ggalluvial, ggdendro, networkD3, htmlwidgets packages)
- Bash shell
- SLURM job scheduler
- ChromHMM and EpiSegMix tools
- BED format segmentation files
- BigWig format signal files

## Output Files

- Classification accuracy summaries (CSV format)
- Feature importance plots (PNG/PDF)
- Model performance comparisons
- Sample-level classification metrics
- Visualization files for manuscript figures
