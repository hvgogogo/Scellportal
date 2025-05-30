# ==============================================================================
# SETUP_SCELLPORTAL.R - Installation and Setup Script
# ==============================================================================

#' Setup ScellPortal Package
#' 
#' This script installs and sets up the ScellPortal package with all dependencies
#' 
#' @author Your Name

# Install required packages if not already installed
install_if_missing <- function(packages) {
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    
    # CRAN packages
    cran_packages <- c("devtools", "dplyr", "ggplot2", "patchwork", "tidyr", 
                       "stringr", "data.table", "ggrepel", "reshape2", "cowplot",
                       "Cairo", "aplot", "ggpubr")
    cran_missing <- intersect(missing_packages, cran_packages)
    if (length(cran_missing) > 0) {
      install.packages(cran_missing)
    }
    
    # Bioconductor packages
    bioc_packages <- c("biomaRt", "clusterProfiler", "org.Hs.eg.db", "msigdbr")
    bioc_missing <- intersect(missing_packages, bioc_packages)
    if (length(bioc_missing) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(bioc_missing)
    }
    
    # GitHub packages
    github_packages <- list(
      "Seurat" = "satijalab/seurat",
      "harmony" = "immunogenomics/harmony",
      "Nebulosa" = "powellgenomicslab/Nebulosa"
    )
    github_missing <- intersect(missing_packages, names(github_packages))
    if (length(github_missing) > 0) {
      for (pkg in github_missing) {
        devtools::install_github(github_packages[[pkg]])
      }
    }
  }
}

# Required packages
required_packages <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "tidyr", "stringr", "data.table",
  "biomaRt", "clusterProfiler", "org.Hs.eg.db", "msigdbr", "harmony", 
  "Nebulosa", "ggrepel", "ggpubr", "reshape2", "cowplot", "Cairo", "aplot"
)

cat("=== ScellPortal Setup ===\n")
cat("Checking and installing required packages...\n")

install_if_missing(required_packages)

cat("✓ All required packages installed\n")
cat("You can now install ScellPortal from GitHub:\n")
cat("devtools::install_github('yourusername/ScellPortal')\n\n")

# ==============================================================================
# EXAMPLE_ANALYSIS.R - Complete Analysis Example
# ==============================================================================

#' ScellPortal Example Analysis
#' 
#' Complete single-cell RNA-seq analysis workflow using ScellPortal
#' Based on the original synovial fibroblast analysis

# Load libraries
library(ScellPortal)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. DATA LOADING AND SETUP
# ==============================================================================

cat("=== ScellPortal Example Analysis ===\n")

# Setup directories
base_path <- "~/scellportal_analysis"
output_dir <- setup_output_dir(file.path(base_path, "results"))

# Load your data (replace with your own data loading)
# For this example, we assume you have a Seurat object called 'obj'
# obj <- readRDS("path/to/your/seurat_object.rds")

# Or create example data for demonstration
if (!exists("obj")) {
  cat("Creating example dataset...\n")
  set.seed(123)
  
  # Create synthetic count matrix
  n_genes <- 2000
  n_cells <- 500
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 3),
    nrow = n_genes,
    dimnames = list(
      paste0("Gene_", 1:n_genes),
      paste0("Cell_", 1:n_cells)
    )
  )
  
  # Create Seurat object
  obj <- CreateSeuratObject(counts = counts)
  
  # Add metadata
  obj$condition <- sample(c("healthy", "disease"), n_cells, replace = TRUE)
  obj$cell_type <- sample(c("SLSFs_CD34pos", "SLSFs_CD34neg", "LLSFs", "SLSFs_IM"), 
                          n_cells, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2))
  obj$cohort <- sample(c("cohort1", "cohort2", "cohort3"), n_cells, replace = TRUE)
  obj$disease_state <- paste(
    sample(c("lean", "obese"), n_cells, replace = TRUE),
    sample(c("healthy", "OA"), n_cells, replace = TRUE),
    sep = "_"
  )
}

cat(sprintf("Loaded dataset with %d genes and %d cells\n", nrow(obj), ncol(obj)))

# ==============================================================================
# 2. PREPROCESSING AND QC
# ==============================================================================

# ==============================================================================
# 7. GENERATE FINAL REPORT
# ==============================================================================

cat("\n=== Generating Final Report ===\n")

# Create analysis summary
analysis_summary <- list(
  dataset_info = list(
    n_genes = nrow(obj),
    n_cells = ncol(obj),
    n_cell_types = length(unique(obj$cell_type)),
    n_conditions = length(unique(obj$condition)),
    conditions = unique(obj$condition),
    cell_types = unique(obj$cell_type)
  ),
  preprocessing = list(
    normalization = "SCTransform",
    integration = "Harmony",
    reduction = "UMAP"
  ),
  differential_expression = if(exists("deg_summary")) {
    list(
      n_comparisons = nrow(deg_summary),
      total_degs = sum(deg_summary$significant_genes),
      cell_types_analyzed = deg_summary$cell_type
    )
  } else NULL,
  files_generated = list(
    plots = list.files(file.path(output_dir, "plots"), recursive = TRUE),
    tables = list.files(file.path(output_dir, "tables"), recursive = TRUE)
  )
)

# Save analysis summary
writeLines(
  jsonlite::toJSON(analysis_summary, pretty = TRUE),
  file.path(output_dir, "analysis_summary.json")
)

# Create README for results
readme_content <- paste0(
  "# ScellPortal Analysis Results\n\n",
  "**Analysis Date:** ", Sys.Date(), "\n",
  "**Dataset:** ", analysis_summary$dataset_info$n_genes, " genes, ", 
  analysis_summary$dataset_info$n_cells, " cells\n",
  "**Cell Types:** ", paste(analysis_summary$dataset_info$cell_types, collapse = ", "), "\n",
  "**Conditions:** ", paste(analysis_summary$dataset_info$conditions, collapse = ", "), "\n\n",
  "## Directory Structure\n\n",
  "- `plots/` - All visualizations\n",
  "  - `umap_overview.pdf` - UMAP plots by condition and cell type\n",
  "  - `violin_key_genes.pdf` - Gene expression violin plots\n",
  "  - `pie_charts.pdf` - Cell composition analysis\n",
  "  - `deg/` - Differential expression plots\n",
  "- `tables/` - Analysis results tables\n",
  "  - `deg/` - Differential expression results\n",
  "  - `enrichment/` - Pathway enrichment results\n",
  "- `qc/` - Quality control plots\n",
  "- `analysis_summary.json` - Complete analysis summary\n\n",
  "## Key Results\n\n"
)

if (exists("deg_summary")) {
  readme_content <- paste0(
    readme_content,
    "### Differential Expression\n",
    "- Analyzed ", nrow(deg_summary), " cell types\n",
    "- Found ", sum(deg_summary$significant_genes), " significant DEGs total\n\n"
  )
}

if (exists("pathway_summary")) {
  readme_content <- paste0(
    readme_content,
    "### Pathway Enrichment\n",
    "- Identified ", pathway_summary$unique_pathways, " significant pathways\n",
    "- See `pathway_summary.txt` for detailed results\n\n"
  )
}

readme_content <- paste0(
  readme_content,
  "## Citation\n\n",
  "If you use these results, please cite ScellPortal:\n",
  "Your Name (2025). ScellPortal: Single-cell RNA-seq Analysis Toolkit. ",
  "https://github.com/yourusername/ScellPortal\n"
)

writeLines(readme_content, file.path(output_dir, "README.md"))

cat("✓ Analysis completed successfully!\n")
cat(sprintf("Results saved to: %s\n", output_dir))
cat("Check README.md for a summary of results\n")

# ==============================================================================
# 8. SESSION INFO AND CLEANUP
# ==============================================================================

# Save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)

cat("\n=== Analysis Summary ===\n")
cat("Dataset:", analysis_summary$dataset_info$n_genes, "genes,", 
    analysis_summary$dataset_info$n_cells, "cells\n")
cat("Cell types:", paste(analysis_summary$dataset_info$cell_types, collapse = ", "), "\n")
cat("Conditions:", paste(analysis_summary$dataset_info$conditions, collapse = ", "), "\n")

if (exists("deg_summary")) {
  cat("DEGs found:", sum(deg_summary$significant_genes), "across", nrow(deg_summary), "cell types\n")
}

cat("Output directory:", output_dir, "\n")
cat("\n✓ ScellPortal analysis pipeline completed!\n")

# ==============================================================================
# CONFIGURATION_TEMPLATE.R - Analysis Configuration Template
# ==============================================================================

#' ScellPortal Analysis Configuration Template
#' 
#' Copy and modify this template for your specific analysis needs

# Analysis configuration
config <- list(
  # Project information
  project_name = "MyProject",
  analysis_date = Sys.Date(),
  author = "Your Name",
  
  # Paths
  data_path = "data/seurat_object.rds",
  output_dir = "results/",
  
  # Preprocessing parameters
  preprocessing = list(
    min_features = 200,
    max_features = 5000,
    mt_cutoff = 20,
    run_sctransform = TRUE,
    integration_method = "harmony",
    batch_var = "orig.ident"
  ),
  
  # Analysis parameters
  analysis = list(
    cluster_resolution = 0.5,
    condition_col = "condition",
    cluster_col = "cell_type",
    test_method = "wilcox",
    logfc_threshold = 0.25,
    min_pct = 0.1
  ),
  
  # Visualization parameters
  visualization = list(
    point_size = 0.5,
    figure_width = 8,
    figure_height = 6,
    add_labels = TRUE
  ),
  
  # Pathway analysis parameters
  pathways = list(
    categories = c("H", "C2"),
    p_cutoff = 0.05,
    min_geneset_size = 10,
    max_geneset_size = 500
  ),
  
  # Color schemes
  colors = list(
    cell_types = c(
      "TypeA" = "#009E73",
      "TypeB" = "#56B4E9",
      "TypeC" = "#E69F00",
      "TypeD" = "#0072B2"
    ),
    conditions = c(
      "control" = "#9DBAD2",
      "treatment" = "#F8BC7E"
    ),
    disease_states = c(
      "healthy" = "#98C897",
      "mild" = "#FE9C9D",
      "severe" = "#CC976B"
    )
  ),
  
  # Genes of interest
  genes_of_interest = list(
    markers = c("CD34", "THY1", "PPARG", "CCL2"),
    inflammatory = c("IL1B", "TNF", "IL6"),
    metabolic = c("PPARG", "ADIPOQ", "LEP")
  )
)

# Save configuration
saveRDS(config, "analysis_config.rds")

# Function to load and validate configuration
load_config <- function(config_path = "analysis_config.rds") {
  if (!file.exists(config_path)) {
    stop("Configuration file not found. Please create one using the template.")
  }
  
  config <- readRDS(config_path)
  
  # Validate required fields
  required_fields <- c("project_name", "data_path", "output_dir")
  missing_fields <- setdiff(required_fields, names(config))
  
  if (length(missing_fields) > 0) {
    stop("Missing required configuration fields: ", paste(missing_fields, collapse = ", "))
  }
  
  return(config)
}

# ==============================================================================
# BATCH_ANALYSIS.R - Batch Processing Script
# ==============================================================================

#' Batch Analysis Script
#' 
#' Process multiple datasets using the same ScellPortal pipeline

run_batch_analysis <- function(config_files, parallel = FALSE) {
  
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    # Parallel processing
    cl <- parallel::makeCluster(min(length(config_files), parallel::detectCores() - 1))
    parallel::clusterEvalQ(cl, library(ScellPortal))
    
    results <- parallel::parLapply(cl, config_files, function(config_file) {
      tryCatch({
        run_single_analysis(config_file)
      }, error = function(e) {
        list(config = config_file, error = e$message, success = FALSE)
      })
    })
    
    parallel::stopCluster(cl)
  } else {
    # Sequential processing
    results <- lapply(config_files, function(config_file) {
      tryCatch({
        run_single_analysis(config_file)
      }, error = function(e) {
        list(config = config_file, error = e$message, success = FALSE)
      })
    })
  }
  
  return(results)
}

run_single_analysis <- function(config_file) {
  cat(sprintf("=== Processing: %s ===\n", config_file))
  
  # Load configuration
  config <- load_config(config_file)
  
  # Load data
  obj <- readRDS(config$data_path)
  
  # Setup output directory
  output_dir <- setup_output_dir(config$output_dir)
  
  # Run analysis pipeline
  # (Include the main analysis steps here, parameterized by config)
  
  cat(sprintf("✓ Completed: %s\n", config$project_name))
  
  return(list(
    config = config_file,
    project = config$project_name,
    output_dir = output_dir,
    success = TRUE
  ))
}

# Example usage:
# config_files <- c("project1_config.rds", "project2_config.rds", "project3_config.rds")
# batch_results <- run_batch_analysis(config_files, parallel = TRUE)

# ==============================================================================
# SCELLPORTAL_WORKFLOW.R - Complete Workflow Function
# ==============================================================================

#' Complete ScellPortal Workflow
#' 
#' Run the entire ScellPortal analysis pipeline with a single function call
#' 
#' @param obj Seurat object
#' @param config Configuration list or path to config file
#' @param run_enrichment Whether to run pathway enrichment analysis
#' @param create_report Whether to generate analysis report
#' @return List with analysis results
#' @export
scellportal_workflow <- function(obj, 
                                config = NULL,
                                run_enrichment = TRUE,
                                create_report = TRUE) {
  
  # Load configuration if path provided
  if (is.character(config)) {
    config <- load_config(config)
  }
  
  # Use default config if none provided
  if (is.null(config)) {
    config <- list(
      output_dir = "scellportal_results/",
      preprocessing = list(min_features = 200, max_features = 5000, mt_cutoff = 20),
      analysis = list(condition_col = "condition", cluster_col = "cell_type"),
      colors = list()
    )
  }
  
  cat("=== ScellPortal Complete Workflow ===\n")
  
  # Setup output directory
  output_dir <- setup_output_dir(config$output_dir)
  
  # 1. Preprocessing
  cat("1. Preprocessing...\n")
  obj <- preprocess_seurat(obj, 
                          min_features = config$preprocessing$min_features %||% 200,
                          max_features = config$preprocessing$max_features %||% 5000,
                          mt_cutoff = config$preprocessing$mt_cutoff %||% 20)
  
  obj <- run_integration(obj, method = "harmony")
  
  # 2. Basic visualizations
  cat("2. Creating visualizations...\n")
  umap_plot <- create_umap_plot(obj, group_by = config$analysis$cluster_col %||% "seurat_clusters")
  ggsave(file.path(output_dir, "plots/umap_overview.pdf"), umap_plot)
  
  # 3. Differential expression
  cat("3. Differential expression analysis...\n")
  if (!is.null(config$analysis$condition_col) && 
      length(unique(obj@meta.data[[config$analysis$condition_col]])) >= 2) {
    
    conditions <- unique(obj@meta.data[[config$analysis$condition_col]])
    deg_results <- find_condition_markers(
      obj,
      condition_col = config$analysis$condition_col,
      group1 = conditions[1],
      group2 = conditions[2]
    )
    
    export_deg_results(deg_results, output_dir = file.path(output_dir, "tables/deg"))
  }
  
  # 4. Pathway enrichment
  if (run_enrichment && exists("deg_results")) {
    cat("4. Pathway enrichment analysis...\n")
    tryCatch({
      enrichment_results <- run_pathway_enrichment(obj, output_dir = file.path(output_dir, "enrichment"))
      export_enrichment_results(enrichment_results, output_dir = file.path(output_dir, "tables/enrichment"))
    }, error = function(e) {
      cat("Pathway enrichment failed:", e$message, "\n")
    })
  }
  
  # 5. Generate report
  if (create_report) {
    cat("5. Generating report...\n")
    # Create analysis summary (code from earlier)
  }
  
  cat("✓ ScellPortal workflow completed!\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(list(
    object = obj,
    output_dir = output_dir,
    config = config
  ))
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

cat("=== ScellPortal Setup Complete ===\n")
cat("Available functions:\n")
cat("- scellportal_workflow(): Complete analysis pipeline\n")
cat("- run_batch_analysis(): Process multiple datasets\n")
cat("- load_config(): Load analysis configuration\n")
cat("\nFor help: ?ScellPortal or vignette('ScellPortal')\n")== Preprocessing and Quality Control ===\n")

# Generate QC plots
qc_plots <- generate_qc_plots(obj, output_dir = file.path(output_dir, "qc"))

# Preprocess the data
obj <- preprocess_seurat(obj, 
                        min_features = 200,
                        max_features = 5000,
                        mt_cutoff = 20,
                        run_sctransform = TRUE)

# Run integration
obj <- run_integration(obj, 
                      method = "harmony",
                      batch_var = "cohort")

cat("✓ Preprocessing completed\n")

# ==============================================================================
# 3. VISUALIZATION
# ==============================================================================

cat("\n=== Creating Visualizations ===\n")

# Define color schemes
cell_type_colors <- c(
  "SLSFs_CD34pos" = "#009E73",
  "SLSFs_CD34neg" = "#56B4E9", 
  "LLSFs" = "#E69F00",
  "SLSFs_IM" = "#0072B2"
)

condition_colors <- c(
  "healthy" = "#9DBAD2",
  "disease" = "#F8BC7E"
)

# Create UMAP plots
umap_celltype <- create_umap_plot(obj, 
                                 group_by = "cell_type",
                                 colors = cell_type_colors,
                                 add_labels = TRUE,
                                 title = "Cell Types")

umap_condition <- create_umap_plot(obj,
                                  group_by = "condition", 
                                  colors = condition_colors,
                                  title = "Condition")

# Combine UMAP plots
umap_combined <- umap_condition / umap_celltype
ggsave(file.path(output_dir, "plots/umap_overview.pdf"), 
       umap_combined, width = 8, height = 10)

# Create violin plots for key genes
key_genes <- c("CD34", "THY1", "PPARG", "CCL2")  # Example genes
violin_plots <- create_violin_plots(obj,
                                   genes = key_genes,
                                   group_by = "cell_type",
                                   colors = cell_type_colors,
                                   add_stats = TRUE)

ggsave(file.path(output_dir, "plots/violin_key_genes.pdf"),
       violin_plots, width = 12, height = 8)

# Create density plots if Nebulosa is available
if (requireNamespace("Nebulosa", quietly = TRUE)) {
  density_plots <- create_density_plots(obj, 
                                       genes = key_genes[1:2],
                                       reduction = "umap")
  
  density_combined <- wrap_plots(density_plots, ncol = 2)
  ggsave(file.path(output_dir, "plots/density_plots.pdf"),
         density_combined, width = 8, height = 4)
}

# Create pie charts for cell composition
pie_charts <- create_pie_charts(obj,
                               group_by = "cell_type",
                               split_by = "condition",
                               colors = cell_type_colors)

pie_combined <- wrap_plots(pie_charts, ncol = 2)
ggsave(file.path(output_dir, "plots/pie_charts.pdf"),
       pie_combined, width = 10, height = 5)

cat("✓ Visualizations created\n")

# ==============================================================================
# 4. DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

cat("\n=== Differential Expression Analysis ===\n")

# Find condition markers
deg_results <- find_condition_markers(
  obj,
  condition_col = "condition",
  group1 = "healthy",
  group2 = "disease",
  cell_types = c("SLSFs_CD34pos", "SLSFs_CD34neg", "LLSFs"),
  cluster_col = "cell_type"
)

# Create DEG plots
deg_plots <- create_deg_plots(deg_results,
                             output_dir = file.path(output_dir, "plots/deg"),
                             plot_types = c("volcano", "waterfall"))

# Export DEG results
export_deg_results(deg_results,
                  output_dir = file.path(output_dir, "tables/deg"))

# Create summary
deg_summary <- summarize_deg_results(deg_results)
write.csv(deg_summary, 
          file.path(output_dir, "tables/deg_summary.csv"), 
          row.names = FALSE)

cat(sprintf("✓ Found DEGs in %d cell types\n", nrow(deg_summary)))

# ==============================================================================
# 5. PATHWAY ENRICHMENT ANALYSIS
# ==============================================================================

cat("\n=== Pathway Enrichment Analysis ===\n")

# Run pathway enrichment (if required packages are available)
if (requireNamespace("clusterProfiler", quietly = TRUE) && 
    requireNamespace("msigdbr", quietly = TRUE)) {
  
  enrichment_results <- run_pathway_enrichment(
    obj,
    conditions = c("healthy", "disease"),
    cluster_col = "cell_type",
    condition_col = "condition",
    categories = c("H", "C2"),
    output_dir = file.path(output_dir, "enrichment")
  )
  
  # Create enrichment dot plot
  if (length(enrichment_results) > 0) {
    enr_dotplot <- create_enrichment_dotplot(
      enrichment_results,
      conditions = c("healthy", "disease"),
      top_pathways = 20
    )
    
    ggsave(file.path(output_dir, "plots/enrichment_dotplot.pdf"),
           enr_dotplot, width = 12, height = 8)
    
    # Export enrichment results
    export_enrichment_results(enrichment_results,
                             output_dir = file.path(output_dir, "tables/enrichment"))
    
    # Create pathway summary report
    pathway_summary <- create_pathway_summary_report(
      enrichment_results,
      output_file = file.path(output_dir, "pathway_summary.txt")
    )
    
    cat("✓ Pathway enrichment analysis completed\n")
  }
} else {
  cat("⚠ Skipping pathway enrichment (missing required packages)\n")
}

# ==============================================================================
# 6. MULTI-CONDITION COMPARISON
# ==============================================================================

cat("\n=== Multi-condition Comparison ===\n")

# If you have multiple disease states
if (length(unique(obj$disease_state)) > 2) {
  multi_results <- multi_condition_comparison(
    obj,
    condition_col = "disease_state",
    cell_types = c("SLSFs_CD34pos", "SLSFs_CD34neg"),
    cluster_col = "cell_type"
  )
  
  # Export multi-condition results
  for (comparison in names(multi_results)) {
    comp_dir <- file.path(output_dir, "tables/multi_condition", comparison)
    export_deg_results(multi_results[[comparison]], output_dir = comp_dir)
  }
  
  cat("✓ Multi-condition comparison completed\n")
}

# ==============================================================================
# 7. GENERATE FINAL REPORT
# ==============================================================================

cat("\n=