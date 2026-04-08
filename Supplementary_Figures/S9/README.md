# Supplementary Figures S9-S12: Specialized Analysis and Model Validation

## Overview
This directory contains scripts and data for **Supplementary Figures S9-S12**, which present specialized analyses including reviewer responses, validation studies, and extended comparisons.

## S9 Contents

### Scripts

#### `pygenome_commands.sh`
Shell script containing pyGenome commands for genomic visualization.
- **Purpose**: Generates landscape plots and genomic feature visualizations
- **Output**: Epilogos plots, chromatin annotation views

## S10 Contents

### Scripts

#### `pygenome_commands.sh`
Shell script with extended pyGenome visualization commands.
- **Purpose**: Alternative visualization approaches for state distributions
- **Features**: Custom annotation tracks, feature overlays

## S11 & S12 Contents

These directories contain additional analysis scripts and extended validation studies supporting the main manuscript findings.

## Workflow Overview

The supplementary figures in S9-S12 provide:
1. **Extended Validations**: Additional model comparison and validation metrics
2. **Alternative Visualizations**: Complement to main figure approaches
3. **Reviewer Responses**: Addressing manuscript review queries
4. **Specialized Analyses**: In-depth investigations of specific features

## Key Analyses

### S9: Extended Model Comparisons
- Multi-tissue validation
- Cross-platform reproducibility
- Extended sample sets

### S10: Visualization Alternatives
- PyGenome epilogos plots
- Landscape-view segmentation
- Feature-centric visualizations

### S11-S12: Validation and Extensions
- Robustness testing
- Parameter sensitivity analysis
- Generalization to new datasets

## Requirements

- PyGenome toolkit
- Bash shell
- Supporting data files from main analysis
- Configuration files for visualization parameters

## Output Files

- PyGenome epilogos plots (PNG)
- Landscape-view visualizations
- Extended validation summaries
- Comparative analysis figures

## Running the Analysis

```bash
# S9-S10: Generate pyGenome visualizations
bash pygenome_commands.sh

# Output files created in designated output directory
```

## Notes

- Exact file listings depend on specific reviewer queries addressed
- See individual figure subdirectories for detailed documentation
- These figures provide supporting evidence for main manuscript claims
- All scripts use sanitized paths compatible with public release
