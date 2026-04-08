# EpiSegMix Manuscript Companion Repository - Nucleic Acids Research

Accompanying scripts and data for **EpiSegMix: Integrating Multiple Epigenomic Signals for Accurate Chromatin State Segmentation**.

## Repository Structure

This repository contains all analysis scripts, configuration files, and key data files used to generate figures for the manuscript published in *Nucleic Acids Research*.

### Main Figures

- **[Figure_1/](./Figure_1/)**: EpiSegMix Classification and Random Forest Analysis
  - Scripts for 154-sample classifier benchmark
  - Random forest model training and evaluation
  - Feature importance analysis
  - Comparison with ChromHMM baseline

- **[Figure_2/](./Figure_2/)**: State Overlap and Heatmap Analysis
  - Concordance analysis between models
  - State-by-state overlap metrics
  - Heatmap visualizations of state patterns
  - ChromHMM vs EpiSegMix comparison

- **[Figure_3/](./Figure_3/)**: Genomic Coverage and Gene Expression Correlation
  - Jaccard index calculations
  - Genome-wide coverage analysis
  - Gene expression correlation with chromatin states
  - Sankey flow diagrams showing state-function relationships

- **[Figure_4/](./Figure_4/)**: Hi-C Chromatin Compartments and Gene Ontology
  - Integration with 3D chromatin structure
  - Hi-C A/B compartment analysis
  - Gene ontology enrichment analysis
  - Functional annotation of chromatin states

### Supplementary Figures

- **[Supplementary_Figures/S1-S3/](./Supplementary_Figures/S1-S3/)**: IHEC Sample Metadata Analysis
  - Cell type distribution analysis
  - ChIP-seq quality control metrics
  - WGBS methylation and coverage statistics

- **[Supplementary_Figures/S4/](./Supplementary_Figures/S4/)**: RNA Integration and Heatmap Variations

- **[Supplementary_Figures/S5/](./Supplementary_Figures/S5/)**: Additional Heatmap Analysis

- **[Supplementary_Figures/S6/](./Supplementary_Figures/S6/)**: Jaccard Indices, Jensen-Shannon Distances, and Coverage Analysis

- **[Supplementary_Figures/S7/](./Supplementary_Figures/S7/)**: HMM Emission Matrix Analysis

- **[Supplementary_Figures/S8/](./Supplementary_Figures/S8/)**: Circular Genome Visualization (Circos Plots)

- **[Supplementary_Figures/S9-S10/](./Supplementary_Figures/S9/)**: PyGenome Extended Visualizations

- **[Supplementary_Figures/S11-S12/](./Supplementary_Figures/S11/)**: Extended Validation and Robustness Testing

- **[Supplementary_Figures/S13/](./Supplementary_Figures/S13/)**: Additional Analyses (formerly "Reviewer Comments")
  - Distribution fitting analysis
  - No-signal state characterization
  - Extended Random Forest feature importance

## Key Features

### Script Types

- **R Scripts (.R, .r)**: Statistical analysis and visualization
  - Data processing and analysis
  - Figure generation
  - Publication-ready output
  
- **Bash Scripts (.sh)**: Workflow orchestration
  - Pipeline automation
  - Parallel job submission (SLURM)
  - Data processing and conversion

### Data Files

- **Configuration Files (.ini, .txt)**: Pipeline parameters
- **CSV/TSV Files**: Analysis results, metadata, classifications
- **Reference Files (Excel)**: State definitions, color schemes, lookup tables

## Getting Started

Each figure directory contains:
1. **README.md**: Detailed documentation of that figure's analysis
2. **Scripts**: Reproducible analysis code (R and/or bash)
3. **Data**: Key input and output files
4. **Configuration**: Parameter files for reproducibility

### Running the Analysis

```bash
# Navigate to specific figure directory
cd Figure_1/

# Read the README for detailed instructions
cat README.md

# Run scripts (adjust paths as needed)
Rscript Random_forest_classifier.R
```

## Requirements

### Software
- **R 3.6+** with packages:
  - tidyverse (dplyr, ggplot2, reshape2)
  - ComplexHeatmap
  - randomForest
  - ggalluvial
  - patchwork
  - circlize
  - gprofiler2

- **Bash shell**
- **SLURM** (optional, for parallel job submission)
- **PyGenome** (optional, for advanced visualizations)

### Data
- IHEC sample metadata
- ChIP-seq signal files (BigWig format)
- WGBS methylation data
- Chromatin state segmentations (BED format)
- Reference genome (hg38)
- Gene annotations (GTF/GFF)
- Hi-C contact data

## File Organization

```
NAR/
├── Figure_1/                    # Classification and feature importance
│   ├── README.md               # Detailed documentation
│   ├── *.R, *.sh              # Analysis scripts
│   └── *.csv, *.txt           # Data and results
├── Figure_2/                    # State overlap analysis
├── Figure_3/                    # Jaccard and expression analysis
├── Figure_4/                    # Gene ontology and Hi-C
├── Supplementary_Figures/       # S1-S13 analyses
│   ├── S1-S3/                 # Metadata analysis
│   ├── S4-S5/                 # Heatmap variations
│   ├── S6-S8/                 # Statistical analysis
│   ├── S9-S10/                # PyGenome visualizations
│   ├── S11-S12/               # Validation studies
│   └── S13/                   # Additional analyses
├── graphical abstract.png       # Manuscript abstract figure
└── README.md                    # This file
```

## Data Paths

All scripts use generic relative paths for reproducibility:
- `./project/` - Project-specific data
- `./data/` - Shared datasets
- `./env/` - Environment/conda paths
- `./resources/` - Reference resources
- `./scratch/` - Temporary working files
- `./output/` - Analysis results

This ensures scripts work across different computing environments without modification.

## Manuscript Abstract

EpiSegMix is a novel chromatin state segmentation approach that integrates multiple epigenomic signals including histone modifications and DNA methylation. By combining hidden Markov models with machine learning validation, EpiSegMix achieves superior accuracy compared to existing methods like ChromHMM across 154 diverse cell types from the International Human Epigenome Consortium (IHEC).

Key findings:
- EpiSegMix and its methylation-aware variant (EpiSegMix+Meth) outperform ChromHMM
- Chromatin states show strong correlation with gene expression
- Integration with Hi-C data reveals 3D-2D chromatin organization coupling
- State-specific gene functions identified through enrichment analysis
- Robust performance across diverse cell types and tissues

## Citation

If you use scripts or analysis approaches from this repository, please cite:

```
[Citation information to be added upon publication]
```

## Related Resources

- **EpiSegMix Tool**: [Repository URL]
- **IHEC Consortium**: http://www.ihec-consortium.org/
- **ChromHMM**: [Reference information]

## Questions and Support

For questions about the manuscript or analysis:
- Review the detailed README.md in each figure directory
- Check script headers for usage information
- Refer to inline script comments for methodological details

## License

[License information to be added]

## Acknowledgments

This work utilized resources from:
- International Human Epigenome Consortium (IHEC)
- DEEP (Deutsche Epigenome Programm)
- ENCODE Project
- [Additional funding/resource acknowledgments]

---

**Last Updated**: 2026-04-08
**Manuscript Status**: In Preparation / Published
**Repository**: Public companion to NAR manuscript
