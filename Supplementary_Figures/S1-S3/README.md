# Supplementary Figures S1-S3: Metadata Analysis and Sample Characterization

## Overview
This directory contains scripts for **Supplementary Figures S1-S3**, which provide comprehensive metadata analysis and characterization of the 154 IHEC samples used in the study, including cell type distributions, ChIP-seq quality control, and methylation coverage statistics.

## Contents

### Scripts

#### `metadata.r`
Comprehensive R script analyzing IHEC sample metadata.
- **Purpose**: Generates metadata visualization figures for supplement
- **Analysis**:
  - Cell type distribution across IHEC samples
  - Sample ontology intermediate level categorization
  - Histone mark availability summary
  - ChIP-seq quality control metrics (total reads, mapping percentage)
  - WGBS methylation and coverage statistics
- **Visualizations**:
  - Bar plots of cell type frequencies
  - Heatmaps of read counts and mapping statistics
  - Coverage plots with error bars
  - Multi-panel figure compositions
  - Color-coded by cell type intermediate classification
- **Output**: PNG and PDF formatted figures suitable for publication

### Data Files

#### Metadata Resources
- **IHEC_metadata_harmonization.v1.2_with_coverage_avg_meth.extended.10x_filtered.sorted.csv**
  - Master metadata file for 154 IHEC samples
  - Contains: Sample IDs, cell types, ontology classifications, coverage metrics
  - Fields: epirr_id, cell_type, harmonized_cell_type, harmonized_sample_ontology_intermediate, mean_methylation, standard_deviation, total_CpGs, CG_coverage

- **epiatlas_chipseq_qc_summary.csv**
  - ChIP-seq quality control metrics
  - Fields: EpiRR ID, antibody (histone mark), Total_Reads, Mapped_Pct, Syn_AUC, Syn_Elbow_Pt, Syn_JSD
  - Covers: H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3, H3K9me3

- **ihec_companion_erirrs.txt**
  - List of samples included in analysis
  - One EpiRR ID per line (without version numbers)

- **IHEC_EpiATLAS_IA_colors_Mar18_2024.json**
  - Color scheme definitions for cell types
  - RGB color codes for consistent visualization
  - Organized by sample ontology intermediate classification
  - Histone mark color assignments

### Workflow

1. **Data Loading**: Load metadata CSV files and JSON color scheme
2. **Data Processing**: 
   - Harmonize sample identifiers
   - Merge metadata from multiple sources
   - Factor categorical variables
   - Calculate summary statistics
3. **Visualization**: Generate cell type and QC plots
4. **Export**: Save publication-quality figures

## Analysis Components

### S1-S3: Cell Type and Metadata Analysis

#### Sample Distribution Analysis
- **Total Samples**: 154 IHEC samples
- **Cell Types**: Categorized by intermediate-level ontology (immune cells, blood, epithelial, fibroblasts, etc.)
- **Sample Frequencies**: Bar charts showing counts per cell type
- **Ontology Mapping**: Shows relationship between cell types and intermediate classifications

#### ChIP-seq Quality Control
- **Histone Marks**: 6 main marks analyzed (H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3, H3K9me3)
- **Quality Metrics**:
  - Total Reads: Total sequencing depth per sample
  - Mapping Percentage: Fraction of reads mapped to reference genome
  - Synthetic AUC: Area under curve for synthetic quality assessment
  - Synthetic Elbow Point: Inflection point in quality curve
  - Jensen-Shannon Divergence: Model quality metric
- **Faceted Plots**: Separate panels for each histone mark
- **Status Indicators**: Color-coded by sample quality tier

#### Methylation and Coverage Statistics
- **WGBS Metrics**:
  - Average Methylation: Mean methylation percentage ± SD across genome
  - Total CpG Sites: Number of covered CpG dinucleotides
  - CpG Coverage: Mean sequencing depth per CpG site (Reads/CpG)
- **Visualization**: Line plots with error bars, coverage distributions
- **Range Checks**: Verifies expected ranges for methylation data

### Color Scheme

Consistent color coding across all figures:
- **By cell type intermediate classification**: Distinct colors for immune, blood, epithelial cells, etc.
- **By histone mark**: Specific color for each of 6 ChIP-seq marks
- **Standardized palette**: Defined in JSON for reproducibility

## Usage Example

```bash
# Generate all S1-S3 figures
Rscript metadata.r

# Output files created:
# - combined_plot.png (combined S1-S3 figure)
# - WGBS_average_methylation_combined.pdf
# - ChIP_QC_combined.pdf
# - cell_type_combined_plot.png
```

## Requirements

- R 3.6+ (with dplyr, stringr, ggplot2, reshape2, scales, jsonlite, patchwork packages)
- Metadata CSV files
- JSON color scheme file
- Writing permissions for output directory

## Output Files

Generated figures:
- **S1-S3 Combined Panel** (PNG): Cell type distribution and sample ontology
- **ChIP-seq QC Panel** (PNG): Histone mark availability and quality metrics
- **WGBS Methylation Panel** (PNG): Methylation and CpG coverage statistics
- **Combined Publication Figure** (PNG/PDF): Multi-panel figure for supplementary material

Each figure includes:
- High resolution (300 DPI)
- Publication-ready formatting
- Labeled axes and legends
- Standardized color scheme
- Sample-level information

## Key Observations

The analysis reveals:
- Diverse cell type representation across 154 IHEC samples
- High-quality ChIP-seq data across all samples and histone marks
- Comprehensive methylation coverage with consistent CpG detection
- Cell type-specific patterns in histone mark distributions
- Quality metrics supporting downstream segmentation analysis

## Data Quality Notes

- **Samples Filtered**: 10x genomics and PCR duplicates removed
- **Sorted**: By sample quality and metadata completeness
- **Validation**: Consistency checks between metadata sources
- **Version Harmonization**: Uses IHEC metadata harmonization v1.2 with extended coverage information

## References

- IHEC (International Human Epigenome Consortium): http://www.ihec-consortium.org/
- Metadata harmonization protocol: Standardized across consortium
- Quality thresholds: Based on IHEC standards for ChIP-seq and WGBS
