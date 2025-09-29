# ScellPortal: Single-cell RNA-seq Analysis Toolkit. 

A comprehensive R package for single-cell RNA-seq data analysis, visualization, and pathway enrichment analysis. (https://hvgogogo.github.io/Scellportal/)

## Features

### 🔬 Core Analysis Functions
- **Data preprocessing and QC**
- **Dimensionality reduction and clustering**
- **Cell type annotation**
- **Differential expression analysis**
- **Pathway enrichment analysis**

### 📊 Visualization Tools
- **UMAP/tSNE plots with annotations**
- **Violin and density plots**
- **Volcano and waterfall plots**
- **Heatmaps and dot plots**
- **Pie charts for cell composition**

### 🧬 Specialized Features
- **Multi-condition comparisons**
- **Gene set enrichment analysis (GSEA)**
- **Custom pathway analysis**
- **Publication-ready figures**

## Package Structure

```
ScellPortal/
├── R/                          # Main R functions
│   ├── preprocessing.R         # Data preprocessing functions
│   ├── visualization.R         # Plotting functions
│   ├── enrichment.R           # Pathway enrichment analysis
│   ├── differential.R         # Differential expression analysis
│   ├── utilities.R            # Helper functions
│   └── themes.R               # Custom ggplot themes
├── man/                       # Documentation files
├── data/                      # Example datasets
├── vignettes/                 # Tutorials and examples
│   ├── basic_workflow.Rmd     # Basic analysis workflow
│   ├── advanced_analysis.Rmd  # Advanced features
│   └── visualization.Rmd      # Plotting examples
├── inst/
│   └── extdata/              # External data files
├── tests/                    # Unit tests
├── DESCRIPTION              # Package metadata
├── NAMESPACE               # Package namespace
└── README.md              # This file
```
