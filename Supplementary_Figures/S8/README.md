# Supplementary Figure S8: Circular Genome Visualization

## Overview
This directory contains scripts for **Supplementary Figure S8**, providing circular (circos) visualizations of chromatin state distributions across the genome.

## Contents

### Scripts

#### `circlize_plot.r`
R script creating circular genome plots and circos diagrams.
- **Purpose**: Visualizes global chromatin state distributions in circular genome format
- **Features**:
  - Circular representation of all chromosomes
  - Radial tracks for state distributions
  - Chord diagrams for state transitions
  - Interactive HTML output with zoom/pan
  - Static PNG/PDF publication format
- **Visualization Types**:
  - Circular state heatmap (all chromosomes)
  - State frequency tracks
  - Transition chord diagrams
  - Model comparison overlays
- **Output**: Publication-ready circos plots

## Analysis Details

### Circos Plot Components

#### Outer Ring: Chromosomal Coordinates
- 22 autosomes + X chromosome
- Megabase pair (Mb) position labels
- Chromosome ideogram representation
- Centromere and p/q arm demarcation

#### Inner Tracks: Chromatin States
- State distributions per chromosome
- Radial distance represents state frequency/presence
- Color-coded by state type:
  - Red: Active chromatin states
  - Orange: Poised/intermediate states
  - Blue: Repressive states
  - Gray: Background/no-signal states
- One track per model compared

#### Chord Connections: State Transitions
- Lines connecting related states
- Line width represents transition frequency
- Visualizes state co-localization
- Identifies chromosomal regions with transitions
- Interactive: click to highlight specific chords

### Biological Patterns Revealed
- **Active states** concentrated in euchromatic regions
- **Repressive states** enriched in heterochromatic regions (pericentromeric)
- **Poised states** at developmental gene loci
- **State transitions** marking regulatory regions
- **Chromosomal variation** in state distributions

## Usage Example

```bash
# Generate circos plots
Rscript circlize_plot.r

# Output files:
# - S8_circos_chromhmm.png
# - S8_circos_episegmix.png
# - S8_circos_episegmix_meth.png
# - S8_circos_interactive.html
# - S8_state_transitions.pdf
```

## Requirements

- R (with circlize, tidyverse, ggplot2 packages)
- BED format segmentation files
- Chromosome size reference file (hg38)
- State definitions and color schemes
- Optional: htmlwidgets for interactive output

## Data Inputs

- Segmentation BED files (one per model)
- Genomic coordinate reference
- Chromosome boundaries
- State color palette definitions

## Output Files

### Static Figures
- **Circos plots (PNG)**: 300 DPI, publication quality
  - One per model (ChromHMM, EpiSegMix, EpiSegMix+Meth)
  - All 24 chromosomes displayed
  - Legend with state definitions
  
- **State transition plots (PDF)**: Chord diagrams
  - Transition frequencies
  - Connection strength indicates co-occurrence

### Interactive Visualization
- **HTML interactive circos**: 
  - Zoom in/out by chromosome
  - Pan across genome
  - Click for detailed information
  - Hover for state information
  - Download as SVG option

### Data Outputs
- **State frequency by chromosome (CSV)**
- **Transition statistics (TSV)**
- **Genome-wide state distribution summary**

## Key Features

### Visual Clarity
- Clean, uncluttered design
- Color-blind friendly palette
- Clear label positioning
- High contrast for black/white printing

### Biological Insight
- Whole-genome perspective (all chromosomes at once)
- Intuitive circular format for genomic data
- Pattern recognition across 3 billion base pairs
- Natural representation of circular genome concept

### Interactive Capabilities
- Explore patterns at different scales
- Link clicking to additional analyses
- Export options for presentations
- Customizable color schemes

## Comparison Across Models

The circos plots enable:
- Visual comparison of ChromHMM vs EpiSegMix
- Assessment of methylation-aware predictions
- Identification of model-specific features
- Validation of state assignments across genome
- Cross-chromosomal pattern consistency

## Interpretation Guidance

**Look for:**
- Consistent state distributions across chromosomes
- Enrichment of specific states on particular chromosomes
- Centromeric vs euchromatic differences
- Sex chromosome (X) specific patterns
- Correlation between models

**Note:**
- Small variations in state positions expected across samples
- Aggregate view masks sample-specific variation
- Rare states may not be visible due to scale
- Interactive version better for exploration
