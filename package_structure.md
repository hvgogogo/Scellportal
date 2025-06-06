# ScellPortal Package Structure

## DESCRIPTION File

```
Package: ScellPortal
Type: Package
Title: Single-cell RNA-seq Analysis Toolkit
Version: 1.0.0
Date: 2025-01-30
Authors@R: c(
    person("Your", "Name", email = "your.email@institution.edu", 
           role = c("aut", "cre"), comment = c(ORCID = "0000-0000-0000-0000"))
    )
Description: A comprehensive toolkit for single-cell RNA-seq data analysis, 
    including preprocessing, visualization, differential expression analysis, 
    and pathway enrichment analysis. Built on the Seurat framework with 
    additional specialized functions for publication-ready visualizations.
License: MIT + file LICENSE
URL: https://github.com/yourusername/ScellPortal
BugReports: https://github.com/yourusername/ScellPortal/issues
Encoding: UTF-8
LazyData: true
RoxygenNote: 7.2.0
Depends: 
    R (>= 4.0.0)
Imports:
    Seurat (>= 4.0.0),
    dplyr,
    ggplot2,
    patchwork,
    tidyr,
    stringr,
    data.table,
    methods,
    stats,
    utils
Suggests:
    biomaRt,
    clusterProfiler,
    org.Hs.eg.db,
    msigdbr,
    harmony,
    Nebulosa,
    ggrepel,
    ggpubr,
    Cairo,
    aplot,
    reshape2,
    cowplot,
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
VignetteBuilder: knitr
Config/testthat/edition: 3
```

## NAMESPACE File

```
# Generated by roxygen2: do not edit by hand

export(annotate_clusters)
export(compare_pathways_across_conditions)
export(convert_gene_ids)
export(create_combined_heatmap)
export(create_density_plots)
export(create_deg_plots)
export(create_dotplot)
export(create_enrichment_dotplot)
export(create_heatmap)
export(create_pathway_abbreviations)
export(create_pathway_summary_report)
export(create_pie_charts)
export(create_umap_plot)
export(create_violin_plots)
export(create_volcano_plot)
export(create_waterfall_plot)
export(density_plot_theme)
export(export_deg_results)
export(export_enrichment_results)
export(filter_enrichment_results)
export(find_all_cluster_markers)
export(find_condition_markers)
export(generate_qc_plots)
export(get_pathway_genesets)
export(get_top_pathways)
export(longitudinal_analysis)
export(multi_condition_comparison)
export(perform_enrichment)
export(pie_chart_theme)
export(preprocess_seurat)
export(prepare_biomart_data)
export(prepare_gene_rankings)
export(process_deg_results)
export(process_enrichment_results)
export(process_pathway_names)
export(run_gsea)
export(run_integration)
export(run_pathway_enrichment)
export(scellportal_theme)
export(setup_output_dir)
export(summarize_deg_results)

import(ggplot2)
import(dplyr)
importFrom(Seurat,DimPlot)
importFrom(Seurat,DotPlot)
importFrom(Seurat,FeatureScatter)
importFrom(Seurat,FetchData)
importFrom(Seurat,FindAllMarkers)
importFrom(Seurat,FindClusters)
importFrom(Seurat,FindMarkers)
importFrom(Seurat,FindNeighbors)
importFrom(Seurat,FindVariableFeatures)
importFrom(Seurat,NormalizeData)
importFrom(Seurat,PercentageFeatureSet)
importFrom(Seurat,PrepSCTFindMarkers)
importFrom(Seurat,RunPCA)
importFrom(Seurat,RunUMAP)
importFrom(Seurat,SCTransform)
importFrom(Seurat,ScaleData)
importFrom(Seurat,VlnPlot)
importFrom(patchwork,wrap_plots)
importFrom(stats,median)
importFrom(utils,combn)
importFrom(utils,head)
importFrom(utils,tail)
```

## LICENSE File

```
MIT License

Copyright (c) 2025 Your Name

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## .Rbuildignore File

```
^.*\.Rproj$
^\.Rproj\.user$
^\.github$
^README\.Rmd$
^LICENSE\.md$
^\.travis\.yml$
^\.gitignore$
^docs$
^_pkgdown\.yml$
^pkgdown$
^\.lintr$
```

## .gitignore File

```
.Rproj.user
.Rhistory
.RData
.Ruserdata
.DS_Store
Thumbs.db
inst/doc
/doc/
/Meta/
*.log
.httr-oauth
*.Rcheck/
*.tar.gz
vignettes/*.html
vignettes/*.pdf
```

## R/zzz.R - Package Startup

```r
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("ScellPortal: Single-cell RNA-seq Analysis Toolkit")
  packageStartupMessage("Version: ", utils::packageDescription("ScellPortal")$Version)
  packageStartupMessage("For help: ?ScellPortal or vignette('ScellPortal')")
  packageStartupMessage("Citation: citation('ScellPortal')")
}

.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op.scellportal <- list(
    scellportal.theme = "default",
    scellportal.colors = "Set2",
    scellportal.verbose = TRUE
  )
  toset <- !(names(op.scellportal) %in% names(op))
  if(any(toset)) options(op.scellportal[toset])
  
  invisible()
}
```

## vignettes/basic_workflow.Rmd

```r
---
title: "Basic ScellPortal Workflow"
author: "Your Name"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Basic ScellPortal Workflow}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)
```

# ScellPortal: Basic Workflow

This vignette demonstrates the basic workflow for single-cell RNA-seq analysis using ScellPortal.

## Installation and Setup

```{r eval=FALSE}
# Install ScellPortal
devtools::install_github("yourusername/ScellPortal")

# Load the package
library(ScellPortal)
library(Seurat)
```

## Loading Data

```{r eval=FALSE}
# Load your Seurat object
# For demonstration, we'll create a simple example
obj <- CreateSeuratObject(counts = matrix(rnorm(1000*100), nrow = 1000))
obj$condition <- sample(c("control", "treatment"), 100, replace = TRUE)
obj$cell_type <- sample(c("TypeA", "TypeB", "TypeC"), 100, replace = TRUE)
```

## Preprocessing

```{r eval=FALSE}
# Basic preprocessing
obj <- preprocess_seurat(obj, 
                        min_features = 200,
                        max_features = 5000,
                        mt_cutoff = 20)

# Integration (if needed)
obj <- run_integration(obj, method = "harmony", batch_var = "orig.ident")
```

## Visualization

```{r eval=FALSE}
# Create UMAP plots
umap_plot <- create_umap_plot(obj, 
                             group_by = "cell_type",
                             add_labels = TRUE,
                             title = "Cell Types")

# Create violin plots for specific genes
violin_plot <- create_violin_plots(obj, 
                                  genes = c("CD34", "THY1"),
                                  group_by = "cell_type")

# Combine plots
library(patchwork)
combined_plot <- umap_plot | violin_plot
```

## Differential Expression Analysis

```{r eval=FALSE}
# Find condition markers
deg_results <- find_condition_markers(
  obj,
  condition_col = "condition",
  group1 = "control",
  group2 = "treatment",
  cell_types = c("TypeA", "TypeB")
)

# Create volcano plots
volcano_plots <- create_deg_plots(deg_results, plot_types = "volcano")
```

## Pathway Enrichment Analysis

```{r eval=FALSE}
# Run pathway enrichment
enrichment_results <- run_pathway_enrichment(
  obj,
  conditions = c("control", "treatment"),
  categories = c("H", "C2")
)

# Create enrichment dot plot
enr_plot <- create_enrichment_dotplot(enrichment_results)
```

## Exporting Results

```{r eval=FALSE}
# Setup output directory
output_dir <- setup_output_dir("results/")

# Export differential expression results
export_deg_results(deg_results, output_dir = file.path(output_dir, "deg/"))

# Export enrichment results
export_enrichment_results(enrichment_results, 
                         output_dir = file.path(output_dir, "enrichment/"))

# Save plots
ggsave(file.path(output_dir, "plots/umap_celltypes.pdf"), 
       umap_plot, width = 8, height = 6)
```

## Session Information

```{r}
sessionInfo()
```
```

## tests/testthat.R

```r
library(testthat)
library(ScellPortal)

test_check("ScellPortal")
```

## tests/testthat/test_preprocessing.R

```r
test_that("preprocess_seurat works correctly", {
  # Create a simple test object
  counts <- matrix(rpois(1000, 10), nrow = 100)
  rownames(counts) <- paste0("Gene_", 1:100)
  colnames(counts) <- paste0("Cell_", 1:10)
  
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj[["percent.mt"]] <- sample(1:30, 10, replace = TRUE)
  
  # Test preprocessing
  result <- preprocess_seurat(obj, min_features = 10, max_features = 50, mt_cutoff = 25)
  
  expect_s4_class(result, "Seurat")
  expect_true("SCT" %in% names(result@assays) || "RNA" %in% names(result@assays))
})

test_that("convert_gene_ids works", {
  skip_if_not_installed("org.Hs.eg.db")
  
  genes <- c("ENSG00000139618", "ENSG00000157764")
  result <- convert_gene_ids(genes, from_type = "ENSEMBL", to_type = "SYMBOL")
  
  expect_is(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_true(nrow(result) >= 0)
})
```

## GitHub Actions Workflow (.github/workflows/R-CMD-check.yaml)

```yaml
on:
  push:
    branches: [main, master]
  pull_request:
    branches: [main, master]

name: R-CMD-check

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}

    name: ${{ matrix.config.os }} (${{ matrix.config.r }})

    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: windows-latest, r: 'release'}
          - {os: macOS-latest, r: 'release'}
          - {os: ubuntu-20.04, r: 'release', rspm: "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"}
          - {os: ubuntu-20.04, r: 'devel', rspm: "https://packagemanager.rstudio.com/cran/__linux__/focal/latest"}

    env:
      R_REMOTES_NO_ERRORS_FROM_WARNINGS: true
      RSPM: ${{ matrix.config.rspm }}
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v2

      - uses: r-lib/actions/setup-r@v1
        with:
          r-version: ${{ matrix.config.r }}

      - uses: r-lib/actions/setup-pandoc@v1

      - name: Query dependencies
        run: |
          install.packages('remotes')
          saveRDS(remotes::dev_package_deps(dependencies = TRUE), ".github/depends.Rds", version = 2)
          writeLines(sprintf("R-%i.%i", getRversion()$major, getRversion()$minor), ".github/R-version")
        shell: Rscript {0}

      - name: Cache R packages
        if: runner.os != 'Windows'
        uses: actions/cache@v2
        with:
          path: ${{ env.R_LIBS_USER }}
          key: ${{ runner.os }}-${{ hashFiles('.github/R-version') }}-1-${{ hashFiles('.github/depends.Rds') }}
          restore-keys: ${{ runner.os }}-${{ hashFiles('.github/R-version') }}-1-

      - name: Install system dependencies
        if: runner.os == 'Linux'
        run: |
          while read -r cmd
          do
            eval sudo $cmd
          done < <(Rscript -e 'writeLines(remotes::system_requirements("ubuntu", "20.04"))')

      - name: Install dependencies
        run: |
          remotes::install_deps(dependencies = TRUE)
          remotes::install_cran("rcmdcheck")
        shell: Rscript {0}

      - name: Check
        env:
          _R_CHECK_CRAN_INCOMING_REMOTE_: false
        run: rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning", check_dir = "check")
        shell: Rscript {0}

      - name: Upload check results
        if: failure()
        uses: actions/upload-artifact@main
        with:
          name: ${{ runner.os }}-r${{ matrix.config.r }}-results
          path: check
```

## inst/CITATION

```
citHeader("To cite ScellPortal in publications use:")

citEntry(
  entry    = "Manual",
  title    = "ScellPortal: Single-cell RNA-seq Analysis Toolkit",
  author   = person("Your Name"),
  year     = "2025",
  note     = "R package version 1.0.0",
  url      = "https://github.com/yourusername/ScellPortal",
  textVersion = paste(
    "Your Name (2025).",
    "ScellPortal: Single-cell RNA-seq Analysis Toolkit.",
    "R package version 1.0.0.",
    "https://github.com/yourusername/ScellPortal"
  )
)
```