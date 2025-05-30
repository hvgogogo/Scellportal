#' ScellPortal: Differential Expression Analysis Functions
#' 
#' Comprehensive differential expression analysis workflows
#' 
#' @author Your Name

# ==============================================================================
# CORE DIFFERENTIAL EXPRESSION FUNCTIONS
# ==============================================================================

#' Find Condition Markers
#' 
#' Find differentially expressed genes between conditions
#' 
#' @param obj Seurat object
#' @param condition_col Column name for conditions
#' @param group1 Vector of conditions for group 1
#' @param group2 Vector of conditions for group 2
#' @param cell_types Vector of cell types to analyze (optional)
#' @param cluster_col Column name for cell types/clusters
#' @param test_method Statistical test method
#' @param logfc_threshold Log fold change threshold
#' @param min_pct Minimum percentage of cells expressing the gene
#' @return List of differential expression results
#' @export
find_condition_markers <- function(obj,
                                  condition_col = "condition",
                                  group1,
                                  group2,
                                  cell_types = NULL,
                                  cluster_col = "cell_type",
                                  test_method = "wilcox",
                                  logfc_threshold = 0.25,
                                  min_pct = 0.1) {
  
  cat("=== Running differential expression analysis ===\n")
  
  # Prepare object
  obj <- Seurat::PrepSCTFindMarkers(obj)
  
  # If no cell types specified, use all
  if (is.null(cell_types)) {
    cell_types <- unique(obj@meta.data[[cluster_col]])
  }
  
  all_results <- list()
  
  for (cell_type in cell_types) {
    cat(sprintf("Analyzing cell type: %s\n", cell_type))
    
    # Subset to current cell type
    obj_subset <- subset(obj, subset = !!sym(cluster_col) == cell_type)
    
    # Create comparison groups
    cells_group1 <- colnames(obj_subset)[obj_subset@meta.data[[condition_col]] %in% group1]
    cells_group2 <- colnames(obj_subset)[obj_subset@meta.data[[condition_col]] %in% group2]
    
    if (length(cells_group1) < 3 || length(cells_group2) < 3) {
      cat(sprintf("  Skipping %s: insufficient cells in one or both groups\n", cell_type))
      next
    }
    
    cat(sprintf("  Group 1 (%s): %d cells\n", paste(group1, collapse = ","), length(cells_group1)))
    cat(sprintf("  Group 2 (%s): %d cells\n", paste(group2, collapse = ","), length(cells_group2)))
    
    # Run differential expression
    tryCatch({
      markers <- Seurat::FindMarkers(
        obj_subset,
        ident.1 = cells_group2,
        ident.2 = cells_group1,
        test.use = test_method,
        logfc.threshold = logfc_threshold,
        min.pct = min_pct,
        verbose = FALSE
      )
      
      # Add gene information
      markers$gene <- rownames(markers)
      if (exists("bm5")) {
        markers$gene_symbol <- bm5$hgnc_symbol[match(markers$gene, bm5$ensembl_gene_id)]
      } else {
        markers$gene_symbol <- markers$gene
      }
      
      # Add comparison information
      markers$cell_type <- cell_type
      markers$comparison <- paste(paste(group2, collapse = "_"), "vs", paste(group1, collapse = "_"))
      
      # Reorder by log fold change
      markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
      
      cat(sprintf("  Found %d significant genes (p_adj < 0.05)\n", 
                  sum(markers$p_val_adj < 0.05)))
      
      all_results[[cell_type]] <- markers
      
    }, error = function(e) {
      cat(sprintf("  Error in %s: %s\n", cell_type, e$message))
    })
  }
  
  cat("✓ Differential expression analysis completed\n\n")
  return(all_results)
}

#' Find All Cluster Markers
#' 
#' Find markers for all clusters in the dataset
#' 
#' @param obj Seurat object
#' @param cluster_col Column name for clusters
#' @param only_pos Only return positive markers
#' @param test_method Statistical test method
#' @param logfc_threshold Log fold change threshold
#' @param min_pct Minimum percentage of cells expressing the gene
#' @return Data frame with all cluster markers
#' @export
find_all_cluster_markers <- function(obj,
                                    cluster_col = "seurat_clusters",
                                    only_pos = TRUE,
                                    test_method = "wilcox",
                                    logfc_threshold = 0.25,
                                    min_pct = 0.25) {
  
  cat("=== Finding markers for all clusters ===\n")
  
  # Prepare object
  obj <- Seurat::PrepSCTFindMarkers(obj)
  
  # Set identities
  Idents(obj) <- cluster_col
  
  # Find all markers
  all_markers <- Seurat::FindAllMarkers(
    obj,
    only.pos = only_pos,
    test.use = test_method,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    verbose = TRUE
  )
  
  # Add gene symbols
  if (exists("bm5")) {
    all_markers$gene_symbol <- bm5$hgnc_symbol[match(all_markers$gene, bm5$ensembl_gene_id)]
  } else {
    all_markers$gene_symbol <- all_markers$gene
  }
  
  cat(sprintf("Found markers for %d clusters\n", length(unique(all_markers$cluster)))
  cat(sprintf("Total markers: %d\n", nrow(all_markers)))
  
  return(all_markers)
}

# ==============================================================================
# ADVANCED DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

#' Multi-Condition Comparison
#' 
#' Perform pairwise comparisons across multiple conditions
#' 
#' @param obj Seurat object
#' @param condition_col Column name for conditions
#' @param conditions Vector of conditions to compare
#' @param cell_types Vector of cell types to analyze
#' @param cluster_col Column name for cell types
#' @return List of comparison results
#' @export
multi_condition_comparison <- function(obj,
                                      condition_col = "condition",
                                      conditions = NULL,
                                      cell_types = NULL,
                                      cluster_col = "cell_type") {
  
  # Get all conditions if not specified
  if (is.null(conditions)) {
    conditions <- unique(obj@meta.data[[condition_col]])
  }
  
  # Get all cell types if not specified
  if (is.null(cell_types)) {
    cell_types <- unique(obj@meta.data[[cluster_col]])
  }
  
  cat(sprintf("=== Multi-condition comparison across %d conditions and %d cell types ===\n",
              length(conditions), length(cell_types)))
  
  # Generate all pairwise comparisons
  comparisons <- combn(conditions, 2, simplify = FALSE)
  
  all_results <- list()
  
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    comp_name <- paste(comp[2], "vs", comp[1])
    
    cat(sprintf("\nComparison %d: %s\n", i, comp_name))
    
    results <- find_condition_markers(
      obj,
      condition_col = condition_col,
      group1 = comp[1],
      group2 = comp[2],
      cell_types = cell_types,
      cluster_col = cluster_col
    )
    
    all_results[[comp_name]] <- results
  }
  
  return(all_results)
}

#' Longitudinal Analysis
#' 
#' Analyze changes across time points or treatment progression
#' 
#' @param obj Seurat object
#' @param time_col Column name for time points
#' @param cell_types Vector of cell types to analyze
#' @param cluster_col Column name for cell types
#' @param reference_time Reference time point for comparison
#' @return List of longitudinal results
#' @export
longitudinal_analysis <- function(obj,
                                 time_col = "timepoint",
                                 cell_types = NULL,
                                 cluster_col = "cell_type",
                                 reference_time = NULL) {
  
  time_points <- sort(unique(obj@meta.data[[time_col]]))
  
  if (is.null(reference_time)) {
    reference_time <- time_points[1]
  }
  
  cat(sprintf("=== Longitudinal analysis with reference: %s ===\n", reference_time))
  
  comparison_times <- setdiff(time_points, reference_time)
  all_results <- list()
  
  for (time_point in comparison_times) {
    comp_name <- paste(time_point, "vs", reference_time)
    cat(sprintf("Analyzing: %s\n", comp_name))
    
    results <- find_condition_markers(
      obj,
      condition_col = time_col,
      group1 = reference_time,
      group2 = time_point,
      cell_types = cell_types,
      cluster_col = cluster_col
    )
    
    all_results[[comp_name]] <- results
  }
  
  return(all_results)
}

# ==============================================================================
# RESULT PROCESSING AND FILTERING
# ==============================================================================

#' Process DEG Results
#' 
#' Process and filter differential expression results
#' 
#' @param markers Data frame with marker results
#' @param fc_cutoff Log2 fold change cutoff
#' @param p_cutoff Adjusted p-value cutoff
#' @param remove_genes Vector of genes to exclude
#' @param top_n Number of top genes to highlight
#' @return List with processed results
#' @export
process_deg_results <- function(markers,
                               fc_cutoff = 1,
                               p_cutoff = 0.05,
                               remove_genes = NULL,
                               top_n = 20) {
  
  # Remove unwanted genes
  if (!is.null(remove_genes)) {
    markers <- markers[!markers$gene %in% remove_genes, ]
  }
  
  # Arrange by p-value
  markers <- markers %>% arrange(p_val_adj)
  
  # Add gene ID column
  markers$gene_id <- rownames(markers)
  
  # Filter significant genes
  up_genes <- markers %>%
    filter(avg_log2FC > fc_cutoff & p_val_adj < p_cutoff)
  
  down_genes <- markers %>%
    filter(avg_log2FC < -fc_cutoff & p_val_adj < p_cutoff)
  
  all_degs <- rbind(up_genes, down_genes)
  
  # Get top genes for labeling
  top_up <- up_genes %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = top_n/2) %>%
    pull(gene_symbol)
  
  top_down <- down_genes %>%
    arrange(avg_log2FC) %>%
    slice_head(n = top_n/2) %>%
    pull(gene_symbol)
  
  features_to_label <- c(top_up, top_down)
  features_to_label <- features_to_label[!is.na(features_to_label)]
  
  result <- list(
    all_genes = markers,
    degs = all_degs,
    up_genes = up_genes,
    down_genes = down_genes,
    features_to_label = features_to_label,
    summary = list(
      total_genes = nrow(markers),
      significant_genes = nrow(all_degs),
      up_regulated = nrow(up_genes),
      down_regulated = nrow(down_genes)
    )
  )
  
  return(result)
}

#' Summarize DEG Results
#' 
#' Create summary statistics for differential expression results
#' 
#' @param deg_results List of DEG results
#' @return Data frame with summary statistics
#' @export
summarize_deg_results <- function(deg_results) {
  
  summary_list <- list()
  
  for (cell_type in names(deg_results)) {
    if (!is.null(deg_results[[cell_type]])) {
      
      markers <- deg_results[[cell_type]]
      
      summary_list[[cell_type]] <- data.frame(
        cell_type = cell_type,
        total_genes = nrow(markers),
        significant_genes = sum(markers$p_val_adj < 0.05),
        up_regulated = sum(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05),
        down_regulated = sum(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05),
        max_logfc = max(abs(markers$avg_log2FC)),
        min_pval = min(markers$p_val_adj),
        stringsAsFactors = FALSE
      )
    }
  }
  
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  return(summary_df)
}

# ==============================================================================
# EXPORT AND VISUALIZATION HELPERS
# ==============================================================================

#' Export DEG Results
#' 
#' Export differential expression results to CSV files
#' 
#' @param deg_results List of DEG results
#' @param output_dir Output directory
#' @param file_prefix Prefix for output files
#' @export
export_deg_results <- function(deg_results,
                               output_dir = "differential_expression/",
                               file_prefix = "DEG") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (cell_type in names(deg_results)) {
    if (!is.null(deg_results[[cell_type]])) {
      
      # Export all results
      all_file <- file.path(output_dir, paste0(file_prefix, "_", cell_type, "_all.csv"))
      write.csv(deg_results[[cell_type]], all_file, row.names = FALSE)
      
      # Export significant results only
      sig_results <- deg_results[[cell_type]] %>%
        filter(p_val_adj < 0.05)
      
      if (nrow(sig_results) > 0) {
        sig_file <- file.path(output_dir, paste0(file_prefix, "_", cell_type, "_significant.csv"))
        write.csv(sig_results, sig_file, row.names = FALSE)
      }
    }
  }
  
  # Create summary file
  summary_df <- summarize_deg_results(deg_results)
  summary_file <- file.path(output_dir, paste0(file_prefix, "_summary.csv"))
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  cat(sprintf("DEG results exported to: %s\n", output_dir))
}

#' Create DEG Plots
#' 
#' Generate comprehensive plots for differential expression results
#' 
#' @param deg_results List of processed DEG results
#' @param output_dir Output directory
#' @param plot_types Vector of plot types to generate
#' @return List of generated plots
#' @export
create_deg_plots <- function(deg_results,
                            output_dir = "deg_plots/",
                            plot_types = c("volcano", "waterfall", "heatmap")) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_plots <- list()
  
  for (cell_type in names(deg_results)) {
    if (!is.null(deg_results[[cell_type]])) {
      
      cat(sprintf("Creating plots for: %s\n", cell_type))
      
      processed_results <- process_deg_results(deg_results[[cell_type]])
      
      cell_plots <- list()
      
      # Volcano plot
      if ("volcano" %in% plot_types) {
        volcano_plot <- create_volcano_plot(
          processed_results$all_genes,
          processed_results$features_to_label
        ) + ggtitle(paste("Volcano Plot -", cell_type))
        
        ggsave(file.path(output_dir, paste0("volcano_", cell_type, ".pdf")),
               volcano_plot, width = 8, height = 6)
        
        cell_plots[["volcano"]] <- volcano_plot
      }
      
      # Waterfall plot
      if ("waterfall" %in% plot_types && nrow(processed_results$degs) > 0) {
        waterfall_plot <- create_waterfall_plot(processed_results$degs) +
          ggtitle(paste("Waterfall Plot -", cell_type))
        
        ggsave(file.path(output_dir, paste0("waterfall_", cell_type, ".pdf")),
               waterfall_plot, width = 10, height = 8)
        
        cell_plots[["waterfall"]] <- waterfall_plot
      }
      
      all_plots[[cell_type]] <- cell_plots
    }
  }
  
  return(all_plots)
}

#' Prepare Gene Rankings
#' 
#' Prepare ranked gene lists for GSEA analysis
#' 
#' @param deg_results DEG results data frame
#' @param ranking_metric Metric to use for ranking ("avg_log2FC", "stat")
#' @param gene_col Column name for gene identifiers
#' @return Named vector of ranked genes
#' @export
prepare_gene_rankings <- function(deg_results,
                                 ranking_metric = "avg_log2FC",
                                 gene_col = "gene_symbol") {
  
  # Remove genes with NA values
  clean_results <- deg_results[!is.na(deg_results[[ranking_metric]]) & 
                               !is.na(deg_results[[gene_col]]), ]
  
  # Create named vector
  gene_rankings <- clean_results[[ranking_metric]]
  names(gene_rankings) <- clean_results[[gene_col]]
  
  # Sort by ranking metric
  gene_rankings <- sort(gene_rankings, decreasing = TRUE)
  
  # Remove duplicates (keep first occurrence)
  gene_rankings <- gene_rankings[!duplicated(names(gene_rankings))]
  
  return(gene_rankings)
}