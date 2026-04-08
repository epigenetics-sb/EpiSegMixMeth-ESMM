# Supplementary Figure S13: Additional Analyses and Extended Studies

## Overview
Directory for **Supplementary Figure S13** (formerly "Reviewers_comments"), containing specialized analyses addressing extended questions and providing additional biological insights beyond main figures.

## Purpose

This directory includes:

1. **Distribution Analysis**:
   - Statistical distributions of chromatin states
   - Fitting analysis for state frequencies
   - Goodness-of-fit testing
   
2. **No-Signal State Analysis**:
   - Characterization of background/no-signal states
   - Bin-level analysis of state transitions
   - Background signal estimation

3. **Feature Importance Analysis**:
   - Extended Random Forest analysis
   - Feature contribution rankings
   - Model interpretation visualizations

4. **Special Case Studies**:
   - Unusual sample analysis
   - Edge case characterization
   - Rare state combinations

## Contents

### Scripts

#### `distribution_fitting.r`
R script analyzing statistical distributions of chromatin states.
- **Purpose**: Fits probability distributions to state frequency data
- **Analysis**:
  - Tests for normality, skewness in state distributions
  - Goodness-of-fit testing (KS, Anderson-Darling)
  - Distribution parameter estimation
  - Visualization of fitted distributions
- **Output**: Distribution plots, goodness-of-fit statistics, parameter tables

#### `No_Signal_state_bins.r`
R script analyzing no-signal/background chromatin states in detail.
- **Purpose**: Characterizes bins assigned to no-signal state
- **Analysis**:
  - Bins with no detectable signal
  - Transition patterns to/from no-signal state
  - Genomic feature associations
  - Sample-specific no-signal prevalence
- **Visualization**: Bar plots, heatmaps of no-signal distribution
- **Output**: No-signal statistics, feature associations

#### `No-singal_state.r`
R script (note: spelling preserved from original) focusing on no-signal state analysis.
- **Purpose**: Extended analysis of background chromatin
- **Features**: State boundary analysis, signal gradient visualization
- **Output**: No-signal state characterization figures

#### `Random_forest_classifier_feature_importance.R`
R script providing extended Random Forest feature importance analysis.
- **Purpose**: Deep dive into feature contributions
- **Analysis**:
  - Variable importance scores for all histone marks
  - Mark-specific discrimination ability
  - Interaction effect analysis
  - Feature selection and ranking
- **Visualization**: 
  - Feature importance bar plots
  - Cumulative importance curves
  - Mark-specific contribution heatmaps
- **Output**: Extended importance rankings, feature contribution matrices

## Key Analyses

### Distribution Analysis
- Tests whether state frequencies follow standard distributions
- Identifies outlier samples or unusual patterns
- Provides statistical context for state prevalence
- Supports claims about model generalization

### No-Signal State Analysis
- **Definition**: Genomic bins with low/no detectable ChIP-seq or WGBS signal
- **Biological Meaning**: Heterochromatin, highly condensed regions
- **Prevalence**: Fraction of genome per sample type
- **Transitions**: How frequently states change to/from no-signal
- **Associations**: Enrichment in repetitive sequences, pericentromeric regions

### Feature Importance Extended Analysis
- **Individual Mark Importance**: Each histone mark's contribution
- **Interaction Terms**: Multi-mark combinations critical for state assignment
- **Redundancy**: Marks that provide similar information
- **Necessity**: Marks essential for model performance

## Output Files

- Distribution fit plots (PNG/PDF)
- No-signal state visualizations (PNG/PDF)
- Feature importance rankings (CSV)
- Extended Random Forest analysis figures
- Statistical test summaries
- Supplementary S13 figure panels
