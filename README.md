# ScellPortal: Single-cell RNA-seq Analysis Toolkit. 

A comprehensive R package for single-cell RNA-seq data analysis, visualization, and pathway enrichment analysis. (https://hvgogogo.github.io/Scellportal/)

## Table of Contents
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Features](#features)
- [Package Structure](#package-structure)
- [Usage Examples](#usage-examples)
- [Documentation](#documentation)
- [Contributing](#contributing)

## Installation

```r
# Install dependencies
source("setup_scellportal.R")

# Install from GitHub
devtools::install_github("hvgogogo/ScellPortal")

# Quick start
library(ScellPortal)
results <- scellportal_workflow(your_seurat_object)
```

### Dependencies

```r
# Required packages
required_packages <- c(
  "Seurat", "biomaRt", "tidyverse", "dplyr", "cowplot",
  "ggplot2", "patchwork", "Cairo", "ggrepel", "Nebulosa",
  "org.Hs.eg.db", "clusterProfiler", "msigdbr", "data.table",
  "reshape2", "ggpubr", "aplot", "stringr"
)

# Install missing packages
install.packages(required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)])
```

## Quick Start

```r
library(ScellPortal)

# Load your Seurat object
data("example_seurat_obj")  # or load your own data

# Basic workflow
obj <- preprocess_seurat(obj)
obj <- run_integration(obj, method = "harmony")
obj <- annotate_clusters(obj)

# Generate comprehensive plots
create_umap_plots(obj, output_dir = "results/")
create_violin_plots(obj, genes = c("CD34", "THY1"), output_dir = "results/")

# Pathway enrichment analysis
enr_results <- run_pathway_enrichment(obj, conditions = c("control", "treatment"))
plot_enrichment_dotplot(enr_results, output_dir = "results/")
```

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

## Usage Examples

### Basic UMAP Visualization

```r
# Create UMAP plots colored by different metadata
p1 <- create_umap_plot(obj, group_by = "cell_type", 
                       colors = custom_colors,
                       title = "Cell Types")

p2 <- create_umap_plot(obj, group_by = "condition",
                       colors = condition_colors,
                       title = "Conditions")

# Combine plots
combined_plot <- p1 / p2
ggsave("umap_overview.pdf", combined_plot, width = 8, height = 10)
```

### Differential Expression Analysis

```r
# Find markers between conditions
deg_results <- find_condition_markers(
  obj, 
  condition_col = "disease_state",
  group1 = c("healthy"),
  group2 = c("disease"),
  cell_types = c("CD34pos_cells", "CD34neg_cells")
)

# Create volcano plots
volcano_plots <- create_volcano_plots(deg_results, 
                                     output_dir = "deg_results/")
```

### Pathway Enrichment

```r
# Run comprehensive pathway analysis
enrichment_results <- run_pathway_enrichment(
  deg_results,
  databases = c("hallmark", "kegg", "reactome"),
  p_cutoff = 0.05
)

# Create enrichment dot plot
enr_dotplot <- plot_enrichment_dotplot(
  enrichment_results,
  top_pathways = 20,
  conditions = c("healthy", "disease")
)

ggsave("pathway_enrichment.pdf", enr_dotplot, width = 12, height = 8)
```

### Custom Visualization Themes

```r
# Apply consistent themes across all plots
my_plots <- list(p1, p2, p3) %>%
  map(~ . + scellportal_theme() + 
      theme(legend.position = "bottom"))
```

## Main Functions Reference

### Data Processing
- `preprocess_seurat()` - Quality control and normalization
- `run_integration()` - Batch correction and integration
- `annotate_clusters()` - Cell type annotation

### Visualization
- `create_umap_plot()` - UMAP visualization
- `create_violin_plots()` - Gene expression violins
- `create_density_plots()` - Gene expression density plots
- `create_heatmap()` - Expression heatmaps
- `create_pie_charts()` - Cell composition analysis

### Analysis
- `find_all_markers()` - Differential expression
- `run_pathway_enrichment()` - Pathway analysis
- `compare_conditions()` - Multi-condition analysis
- `run_gsea()` - Gene set enrichment analysis

### Utilities
- `convert_gene_ids()` - Gene ID conversion
- `prepare_pathway_data()` - Pathway database preparation
- `export_results()` - Export analysis results

## Configuration

Create a configuration file for your analysis:

```r
# config.R
config <- list(
  # Paths
  base_path = "~/data/scrnaseq/",
  output_dir = "results/",
  
  # Analysis parameters
  resolution = 0.5,
  min_cells = 3,
  min_features = 200,
  
  # Visualization
  point_size = 0.5,
  figure_width = 8,
  figure_height = 6,
  
  # Colors
  cell_type_colors = c("#9C6BA3", "#FE9C9D", "#98C897"),
  condition_colors = c("#9DBAD2", "#F8BC7E", "#CC976B")
)
```

## Documentation

Detailed documentation is available in the `vignettes/` directory:

- **Basic Workflow**: Step-by-step analysis guide
- **Advanced Analysis**: Complex multi-condition studies
- **Visualization Guide**: Creating publication-ready figures
- **Pathway Analysis**: Comprehensive enrichment analysis

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup

```bash
git clone https://github.com/hvgogogo/ScellPortal.git
cd ScellPortal
R -e "devtools::load_all(); devtools::check()"
```

## Citation

If you use ScellPortal in your research, please cite:

```
@software{scellportal2025,
  title = {ScellPortal: Single-cell RNA-seq Analysis Toolkit},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/hvgogogo/ScellPortal}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Your Name
- **Email**: your.email@institution.edu
- **GitHub**: [@hvgogogo](https://github.com/hvgogogo)

## Acknowledgments

- Built on the excellent [Seurat](https://satijalab.org/seurat/) framework
- Pathway analysis powered by [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- Visualization enhanced with [ggplot2](https://ggplot2.tidyverse.org/) and [patchwork](https://patchwork.data-imaginist.com/)
