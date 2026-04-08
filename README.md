

# Integrated flexible DNA methylation-chromatin segmentation modeling enhances epigenomic state annotation

![Graphical Abstract](graphical%20abstract.png)

**Manuscript:** Available on bioRxiv  
**DOI:** https://doi.org/10.1101/2025.07.25.666820

This repository contains analysis scripts and supplementary data supporting the manuscript, organized by figure.

If you want to use EpiSegMix or EpiSegMixMeth for your own segmentation analysis, please refer to https://gitlab.com/rahmannlab/episegmix.

---

## Citation

If you use the scripts from this repository:

```bibtex
@article{,
  title={Integrated flexible DNA methylation-chromatin segmentation modeling enhances epigenomic state annotation},
  journal={bioRxiv},
  doi={10.1101/2025.07.25.666820},
  year={2025}
}
```

If you use the EpiSegMix tool, please also cite the original publication:

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
├── Figure_1/                  # Classification and Random Forest Analysis
├── Figure_2/                  # State Overlap and Heatmap Analysis
├── Figure_3/                  # Genomic Coverage and Gene Expression Analysis
├── Figure_4/                  # Hi-C Chromatin Compartments and Gene Ontology
├── Supplementary_Figures/     # S1-S13 Supplementary Analyses
├── graphical abstract.png
└── README.md
```

---

## Notes

Scripts are written in R and Bash and use relative paths (e.g., `./project/`, `./data/`). Full reproducibility requires the complete IHEC dataset, which is not included here due to size. Scripts are provided for methodological transparency and reference.