# EpiSegMix Manuscript Companion Repository

## Integrated flexible DNA methylation-chromatin segmentation modeling enhances epigenomic state annotation

![Graphical Abstract](graphical%20abstract.png)

**Manuscript:** Currently available on bioRxiv  
**DOI:** https://doi.org/10.1101/2025.07.25.666820

This repository contains analysis scripts and supplementary data supporting the manuscript. All computational analyses are organized by figure with complete code for reproducibility.

---

## Citation

### If you use the scripts from this repository:

```bibtex
@article{,
  title={Integrated flexible DNA methylation-chromatin segmentation modeling enhances epigenomic state annotation},
  journal={bioRxiv},
  doi={10.1101/2025.07.25.666820},
  year={2025}
}
```

### If you use the EpiSegMix tool:

Please cite both the manuscript above AND the original EpiSegMix publication:

```bibtex
@article{EpiSegMix,
  title={EpiSegMix: A flexible mixture model approach for epigenomic segmentation},
  journal={Oxford Bioinformatics},
  year={2024}
}
```

---

## Repository Structure

```
NAR/
├── Figure_1/                        # Classification and Random Forest Analysis
├── Figure_2/                        # State Overlap and Heatmap Analysis  
├── Figure_3/                        # Genomic Coverage and Gene Expression Analysis
├── Figure_4/                        # Hi-C Chromatin Compartments and Gene Ontology
├── Supplementary_Figures/           # S1-S13 Supplementary Analyses
│   ├── S1-S3/                      # IHEC Sample Metadata Analysis
│   ├── S4-S12/                     # Additional Analyses
│   └── S13_Additional_Analyses/    # Extended Validation
├── graphical abstract.png           # Manuscript graphical abstract
└── README.md                        # This file
```

---

## Scripts and Analysis Code

This repository contains **28 analysis scripts** across R and Bash:

- **R Scripts**: Statistical analysis and visualization
- **Bash Scripts**: Workflow orchestration and pipeline automation
- **Configuration Files**: Pipeline parameters and settings
- **Data Files**: Analysis results and metadata

All scripts use generic relative paths (e.g., `./project/`, `./data/`, `./resources/`) for cross-platform compatibility.

---

## Requirements

### Software
- **R 3.6+** with core packages (tidyverse, ggplot2, ComplexHeatmap, randomForest)
- **Bash shell**
- **SLURM** (optional, for parallel job submission)

### Data
Full reproducibility requires the complete IHEC dataset (large size). Scripts are provided for reference and methodological transparency.

---

## Key Features

- **Comprehensive Analysis**: 28 analysis scripts covering all figures and supplementary analyses
- **Flexible Paths**: All scripts use sanitized, generic paths for reproducibility
- **Complete Documentation**: Detailed comments and documentation in all scripts
- **Validated Methods**: Code integrity certified—only paths sanitized, all logic preserved

---

## Methodology Overview

This manuscript presents EpiSegMix, a flexible chromatin state segmentation approach that integrates multiple epigenomic signals (histone modifications, DNA methylation) using hidden Markov models validated with machine learning. The analysis demonstrates superior performance compared to existing methods across 154 diverse cell types from the International Human Epigenome Consortium (IHEC).

**Key Analyses:**
1. Classifier benchmark and feature importance
2. State concordance and overlap analysis
3. Genomic coverage and gene expression correlation
4. 3D-2D chromatin organization integration
5. Gene ontology enrichment and functional annotation

---

## Citation and Support

For citation formats and additional information, please refer to the bioRxiv preprint:  
https://doi.org/10.1101/2025.07.25.666820

For the EpiSegMix tool and method details, consult the Oxford Bioinformatics publication.

---

**Repository Status**: Production-ready for GitHub release  
**Manuscript Status**: bioRxiv preprint  
**Last Updated**: 2026-04-08
