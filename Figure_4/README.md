# Figure 4: Hi-C Compartments and Gene Ontology Analysis

## Overview
This directory contains scripts for **Figure 4**, which integrates Hi-C chromatin compartment data with segmentation results and performs gene ontology enrichment analysis on state-specific genes.

## Contents

### Scripts

#### `Hi-C_overlap_compartments.R`
R script analyzing overlap between chromatin states and Hi-C compartments.
- **Purpose**: Correlates EpiSegMix/ChromHMM states with 3D chromatin organization
- **Analysis**:
  - Maps segmentation states to A/B compartments from Hi-C data
  - Calculates compartment enrichment per state
  - Identifies state-specific 3D chromatin organization patterns
  - Generates enrichment heatmaps
- **Visualization**: Compartment × state heatmaps, enrichment plots
- **Output**: Correlation matrices, compartment assignment predictions

#### `gene_ontology_gprofiler.r`
Comprehensive R script for gene ontology (GO) enrichment analysis.
- **Purpose**: Identifies biological functions enriched in genes within each chromatin state
- **Analysis**:
  - Groups genes by chromatin state
  - Performs GO enrichment using g:Profiler database
  - Tests biological process, molecular function, and cellular component ontologies
  - Calculates enrichment p-values and effect sizes
- **Visualization**: GO term bubble plots, term enrichment rankings
- **Output**: Enriched GO terms, ontology heatmaps, interactive results

#### `sankey_overlap.r`
R script creating Sankey diagrams showing relationships between states and biological functions.
- **Purpose**: Visualizes connections between chromatin states and their functional annotations
- **Features**:
  - Shows flow from chromatin states to enriched biological processes
  - Node sizing based on gene counts
  - Color-coded by functional category
  - Interactive HTML widget format
- **Output**: Interactive Sankey visualizations, static figures

#### `pygenome_commands.sh`
Shell script with commands for pyGenome-based genomic feature visualization.
- **Purpose**: Generates landscape plots of chromatin states and functional regions
- **Usage**: Command reference for genomic annotation and feature visualization

### Workflow

1. **Hi-C Integration**: Run `Hi-C_overlap_compartments.R` to link states with 3D structure
2. **GO Enrichment**: Execute `gene_ontology_gprofiler.r` to identify biological processes
3. **Functional Visualization**: Use `sankey_overlap.r` for integrated state-function mapping
4. **Figure Generation**: Create publication-quality panels from analysis outputs

## Analysis Details

### Hi-C Compartment Analysis
- **Data Source**: Hi-C contact matrices at 10kb resolution
- **Compartments**: 
  - **A compartments**: Euchromatic, gene-rich regions
  - **B compartments**: Heterochromatic, gene-poor regions
- **State Enrichment**: Calculates fraction of each state in A vs B compartments
- **3D-2D Integration**: Links local chromatin states to global 3D organization

### Gene Ontology Enrichment
- **Resource**: g:Profiler database (comprehensive gene ontology coverage)
- **Ontologies Tested**:
  - Gene Ontology (GO) Biological Process
  - Gene Ontology (GO) Molecular Function
  - Gene Ontology (GO) Cellular Component
  - KEGG pathways
  - Reactome pathways
- **Statistical Testing**: Hypergeometric test with multiple testing correction
- **Significance**: FDR-adjusted p-value < 0.05

### Sankey Visualization
- **Node 1**: Chromatin states (11 categories)
- **Node 2**: Top enriched GO terms or pathways
- **Flow Thickness**: Number of genes connecting state to function
- **Color Scheme**: Distinct colors for biological process categories

## Usage Example

```bash
# Analyze Hi-C compartment overlap
Rscript Hi-C_overlap_compartments.R

# Perform GO enrichment analysis
Rscript gene_ontology_gprofiler.r

# Create functional relationship visualizations
Rscript sankey_overlap.r

# Optional: Generate pyGenome landscape plots
bash pygenome_commands.sh
```

## Requirements

- R (with tidyverse, ggplot2, ggalluvial, igraph packages)
- g:Profiler API access (installed via gprofiler2 R package)
- Hi-C contact matrix file (HDF5 or text format)
- Gene annotation with coordinates (GTF/GFF)
- Chromatin state assignments from Figure 1

## Data Inputs

- Segmentation BED files (ChromHMM, EpiSegMix, EpiSegMix+Meth)
- Hi-C compartment assignments (per genomic region)
- Gene expression or detection information
- GO database annotations
- Genomic coordinate reference files

## Output Files

- State × compartment enrichment matrices (CSV)
- Compartment overlap plots (PNG/PDF)
- GO enrichment results (TSV/CSV format)
- GO term bubble plots (PNG/PDF)
- Sankey diagrams (HTML interactive, PNG static)
- Functional annotation heatmaps
- Manuscript-ready figure panels

## Key Findings

The analysis reveals:
- Strong enrichment of specific chromatin states in A/B compartments
- State-specific gene functions related to biological processes
- Distinct 3D chromatin organization of different segmentation states
- Enriched pathways in state-specific gene sets
- Correlation between chromatin state, 3D structure, and biological function
- Cell type-specific functional state assignments

## Interpretation Notes

- **Active states** typically show enrichment in:
  - A compartments
  - GO terms related to transcription, translation
  - Developmental and cell differentiation processes

- **Repressive states** typically show enrichment in:
  - B compartments  
  - GO terms related to transcriptional silencing
  - Developmental repression and differentiation pathways

- **Poised states** show intermediate characteristics
