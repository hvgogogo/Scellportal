#' ScellPortal: Pathway Enrichment Analysis Functions
#' 
#' Comprehensive pathway and gene set enrichment analysis
#' 
#' @author Your Name

# ==============================================================================
# PATHWAY DATABASE PREPARATION
# ==============================================================================

#' Get Pathway Gene Sets
#' 
#' Retrieve pathway gene sets from MSigDB
#' 
#' @param gene_type Type of gene IDs ("SYMBOL" or "ENTREZID")
#' @param species Species name
#' @param categories Vector of MSigDB categories to include
#' @return List of pathway gene sets
#' @export
get_pathway_genesets <- function(gene_type = "SYMBOL", 
                                species = "Homo sapiens",
                                categories = c("H", "C2", "C5")) {
  
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr package required for pathway analysis")
  }
  
  cat("=== Retrieving pathway gene sets ===\n")
  
  gene_col <- ifelse(gene_type == "ENTREZID", "entrez_gene", "gene_symbol")
  pathway_sets <- list()
  
  # Hallmark pathways (H)
  if ("H" %in% categories) {
    cat("Retrieving Hallmark pathways...\n")
    hallmark_df <- msigdbr::msigdbr(species = species, category = "H")
    pathway_sets[["hallmark"]] <- hallmark_df %>%
      dplyr::select(gs_name, !!sym(gene_col))
  }
  
  # C2 pathways (curated gene sets)
  if ("C2" %in% categories) {
    cat("Retrieving C2 pathways...\n")
    c2_df <- msigdbr::msigdbr(species = species, category = "C2")
    
    # KEGG pathways
    pathway_sets[["kegg"]] <- c2_df %>%
      dplyr::filter(gs_subcat == "CP:KEGG") %>%
      dplyr::select(gs_name, !!sym(gene_col))
    
    # REACTOME pathways
    pathway_sets[["reactome"]] <- c2_df %>%
      dplyr::filter(gs_subcat == "CP:REACTOME") %>%
      dplyr::select(gs_name, !!sym(gene_col))
    
    # WikiPathways
    pathway_sets[["wiki"]] <- c2_df %>%
      dplyr::filter(gs_subcat == "CP:WIKIPATHWAYS") %>%
      dplyr::select(gs_name, !!sym(gene_col))
  }
  
  # C5 Gene Ontology pathways
  if ("C5" %in% categories) {
    cat("Retrieving C5 Gene Ontology pathways...\n")
    c5_df <- msigdbr::msigdbr(species = species, category = "C5")
    
    # GO Biological Process
    pathway_sets[["go_bp"]] <- c5_df %>%
      dplyr::filter(gs_subcat == "GO:BP") %>%
      dplyr::select(gs_name, !!sym(gene_col))
    
    # GO Cellular Component
    pathway_sets[["go_cc"]] <- c5_df %>%
      dplyr::filter(gs_subcat == "GO:CC") %>%
      dplyr::select(gs_name, !!sym(gene_col))
    
    # GO Molecular Function
    pathway_sets[["go_mf"]] <- c5_df %>%
      dplyr::filter(gs_subcat == "GO:MF") %>%
      dplyr::select(gs_name, !!sym(gene_col))
  }
  
  # Create summary
  pathway_summary <- data.frame(
    pathway_type = names(pathway_sets),
    gene_sets_count = sapply(pathway_sets, function(x) length(unique(x$gs_name))),
    total_genes = sapply(pathway_sets, nrow)
  )
  
  cat("Pathway sets summary:\n")
  print(pathway_summary)
  
  result <- list(
    pathway_sets = pathway_sets,
    summary = pathway_summary,
    gene_type_used = gene_type,
    species_used = species
  )
  
  cat("✓ Pathway gene sets retrieval completed\n\n")
  return(result)
}

# ==============================================================================
# ENRICHMENT ANALYSIS
# ==============================================================================

#' Perform Enrichment Analysis
#' 
#' Run over-representation analysis using clusterProfiler
#' 
#' @param genes_df Data frame with converted gene IDs
#' @param pathway_sets List of pathway gene sets
#' @param pvalue_cutoff P-value cutoff for significance
#' @param min_geneset_size Minimum gene set size
#' @param max_geneset_size Maximum gene set size
#' @return List of enrichment results
#' @export
perform_enrichment <- function(genes_df, 
                              pathway_sets, 
                              pvalue_cutoff = 0.05,
                              min_geneset_size = 10,
                              max_geneset_size = 500) {
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package required for enrichment analysis")
  }
  
  cat("=== Performing enrichment analysis ===\n")
  
  enr_results <- list()
  
  for (pathway_name in names(pathway_sets)) {
    cat(sprintf("Processing pathway set: %s\n", pathway_name))
    
    # Run enrichment
    enr_results[[pathway_name]] <- clusterProfiler::enricher(
      genes_df$converted_id,
      TERM2GENE = pathway_sets[[pathway_name]],
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = "BH",
      minGSSize = min_geneset_size,
      maxGSSize = max_geneset_size
    )
    
    # Add additional metrics
    if (nrow(enr_results[[pathway_name]]@result) > 0) {
      result_df <- enr_results[[pathway_name]]@result
      
      # Extract gene set size from BgRatio
      result_df$genesets_num <- sapply(result_df$BgRatio, function(x) {
        as.numeric(unlist(strsplit(x, "/"))[1])
      })
      
      # Calculate gene ratio
      result_df$GeneRatio <- result_df$Count / result_df$genesets_num
      
      # Update the result
      enr_results[[pathway_name]]@result <- result_df
      
      significant_count <- nrow(result_df)
    } else {
      significant_count <- 0
    }
    
    cat(sprintf("  - Significant terms found: %d\n", significant_count))
  }
  
  cat("✓ Enrichment analysis completed\n\n")
  return(enr_results)
}

#' Process Enrichment Results
#' 
#' Process and format enrichment analysis results
#' 
#' @param enr_results List of enrichment results
#' @param cluster_id Cluster identifier
#' @param condition Condition name
#' @param p_cutoff P-value cutoff for significance
#' @return List with processed results
#' @export
process_enrichment_results <- function(enr_results, 
                                      cluster_id, 
                                      condition,
                                      p_cutoff = 0.01) {
  
  cat(sprintf("=== Processing results for %s - %s ===\n", condition, cluster_id))
  
  # Combine all pathway results
  all_paths <- data.frame()
  
  for (pathway_name in names(enr_results)) {
    if (nrow(enr_results[[pathway_name]]@result) > 0) {
      pathway_results <- enr_results[[pathway_name]]@result
      pathway_results$pathway_source <- pathway_name
      all_paths <- rbind(all_paths, pathway_results)
    }
  }
  
  if (nrow(all_paths) > 0) {
    # Process pathway names
    all_paths <- process_pathway_names(all_paths)
    
    # Filter significant pathways
    significant_paths <- all_paths %>%
      filter(p.adjust <= p_cutoff) %>%
      arrange(p.adjust)
    
    cat(sprintf("Total pathways tested: %d\n", nrow(all_paths)))
    cat(sprintf("Significant pathways (p.adj <= %g): %d\n", p_cutoff, nrow(significant_paths)))
    
    # Get top pathways
    if (nrow(significant_paths) > 0) {
      top_pathways <- head(significant_paths$pathway, 10)
      cat("Top 10 significant pathways:\n")
      for (i in 1:length(top_pathways)) {
        cat(sprintf("  %d. %s\n", i, top_pathways[i]))
      }
    }
    
    result <- list(
      significant_paths = significant_paths,
      all_paths = all_paths
    )
    
    cat("✓ Results processing completed\n\n")
    return(result)
  } else {
    cat("No enrichment results found for this cluster\n\n")
    return(NULL)
  }
}

#' Process Pathway Names
#' 
#' Clean and format pathway names for visualization
#' 
#' @param df Data frame with pathway results
#' @return Processed data frame
#' @export
process_pathway_names <- function(df) {
  
  cat("Processing pathway names and descriptions...\n")
  
  original_rows <- nrow(df)
  data.table::setDT(df)
  
  # Process pathway names
  df$pathway <- df$Description
  df[, pathway := gsub("_", " ", pathway)]
  df$pathway <- sub(" ", ":", df$pathway, fixed = TRUE)
  df$type <- df$pathway
  df$pathway <- sub("^[^:]*:", "", df$pathway)
  
  # Clean pathway names
  df <- df %>%
    mutate(
      pathway = tolower(pathway),
      pathway = stringr::str_replace(pathway, "^(\\w)", function(x) toupper(x))
    )
  
  df$type <- sub(":.*$", "", df$type)
  
  # Reorder columns
  df <- df %>%
    dplyr::select(
      type, pathway, p.adjust, geneID, Count,
      GeneRatio, genesets_num, everything()
    )
  
  cat(sprintf("Processed %d pathway names\n", original_rows))
  
  return(df)
}

# ==============================================================================
# COMPREHENSIVE WORKFLOW FUNCTIONS
# ==============================================================================

#' Run Pathway Enrichment Workflow
#' 
#' Complete pathway enrichment analysis workflow
#' 
#' @param obj Seurat object
#' @param conditions Vector of conditions to analyze
#' @param cluster_col Column name for clusters
#' @param condition_col Column name for conditions
#' @param categories Pathway categories to include
#' @param output_dir Output directory for results
#' @return List of enrichment results
#' @export
run_pathway_enrichment <- function(obj,
                                  conditions = NULL,
                                  cluster_col = "seurat_clusters",
                                  condition_col = "condition",
                                  categories = c("H", "C2"),
                                  output_dir = "enrichment_results/") {
  
  # Setup output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get pathway databases
  pathway_data <- get_pathway_genesets(categories = categories)
  
  # If no conditions specified, use all unique conditions
  if (is.null(conditions)) {
    conditions <- unique(obj@meta.data[[condition_col]])
  }
  
  all_results <- list()
  
  for (condition in conditions) {
    cat(sprintf("\n=== Processing condition: %s ===\n", condition))
    
    # Subset object for current condition
    obj_subset <- subset(obj, subset = !!sym(condition_col) == condition)
    obj_subset <- Seurat::PrepSCTFindMarkers(obj_subset)
    
    # Set cluster identities
    Idents(obj_subset) <- cluster_col
    
    # Find markers for each cluster
    cat("Finding cluster markers...\n")
    all_markers <- Seurat::FindAllMarkers(
      obj_subset,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    # Add gene symbols if using Ensembl IDs
    if (exists("bm5")) {
      all_markers$gene_symbol <- bm5$hgnc_symbol[match(all_markers$gene, bm5$ensembl_gene_id)]
    } else {
      all_markers$gene_symbol <- all_markers$gene
    }
    
    condition_results <- list()
    
    # Process each cluster
    clusters <- unique(all_markers$cluster)
    for (cluster in clusters) {
      cat(sprintf("Processing cluster: %s\n", cluster))
      
      # Get significant genes for this cluster
      cluster_genes <- all_markers %>%
        filter(cluster == !!cluster, p_val_adj < 0.05) %>%
        pull(gene)
      
      if (length(cluster_genes) > 0) {
        # Convert gene IDs
        gene_conversion <- convert_gene_ids(cluster_genes)
        
        if (nrow(gene_conversion$converted_genes) > 0) {
          # Run enrichment analysis
          enrichment_results <- perform_enrichment(
            gene_conversion$converted_genes,
            pathway_data$pathway_sets
          )
          
          # Process results
          processed_results <- process_enrichment_results(
            enrichment_results,
            cluster,
            condition
          )
          
          condition_results[[cluster]] <- processed_results
        }
      }
    }
    
    all_results[[condition]] <- condition_results
    
    # Save condition results
    save_file <- file.path(output_dir, paste0(condition, "_enrichment.RData"))
    save(condition_results, file = save_file)
    cat(sprintf("Results saved to: %s\n", save_file))
  }
  
  # Save combined results
  save(all_results, file = file.path(output_dir, "all_enrichment_results.RData"))
  
  return(all_results)
}

#' Run GSEA Analysis
#' 
#' Gene Set Enrichment Analysis using ranked gene lists
#' 
#' @param ranked_genes Named vector of genes with ranking metric (e.g., log2FC)
#' @param pathway_sets List of pathway gene sets
#' @param nperm Number of permutations
#' @param min_geneset_size Minimum gene set size
#' @param max_geneset_size Maximum gene set size
#' @return GSEA results object
#' @export
run_gsea <- function(ranked_genes,
                    pathway_sets,
                    nperm = 1000,
                    min_geneset_size = 15,
                    max_geneset_size = 500) {
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package required for GSEA")
  }
  
  cat("=== Running GSEA analysis ===\n")
  
  # Sort genes by ranking metric
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  gsea_results <- list()
  
  for (pathway_name in names(pathway_sets)) {
    cat(sprintf("Running GSEA for: %s\n", pathway_name))
    
    gsea_results[[pathway_name]] <- clusterProfiler::GSEA(
      geneList = ranked_genes,
      TERM2GENE = pathway_sets[[pathway_name]],
      nPerm = nperm,
      minGSSize = min_geneset_size,
      maxGSSize = max_geneset_size,
      pvalueCutoff = 1.0,  # Keep all results for filtering later
      verbose = FALSE
    )
  }
  
  cat("✓ GSEA analysis completed\n")
  return(gsea_results)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Export Enrichment Results
#' 
#' Export enrichment results to CSV files
#' 
#' @param enrichment_results List of enrichment results
#' @param output_dir Output directory
#' @param file_prefix Prefix for output files
#' @export
export_enrichment_results <- function(enrichment_results,
                                     output_dir = "results/",
                                     file_prefix = "enrichment") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (condition in names(enrichment_results)) {
    for (cluster in names(enrichment_results[[condition]])) {
      if (!is.null(enrichment_results[[condition]][[cluster]])) {
        
        # Export significant pathways
        sig_file <- file.path(output_dir, 
                             paste0(file_prefix, "_", condition, "_", cluster, "_significant.csv"))
        
        write.csv(enrichment_results[[condition]][[cluster]]$significant_paths,
                 sig_file, row.names = FALSE)
        
        # Export all pathways
        all_file <- file.path(output_dir, 
                             paste0(file_prefix, "_", condition, "_", cluster, "_all.csv"))
        
        write.csv(enrichment_results[[condition]][[cluster]]$all_paths,
                 all_file, row.names = FALSE)
      }
    }
  }
  
  cat(sprintf("Enrichment results exported to: %s\n", output_dir))
}

#' Create Pathway Abbreviations
#' 
#' Create abbreviated pathway names for visualization
#' 
#' @param pathway_names Vector of pathway names
#' @return Named vector of abbreviations
#' @export
create_pathway_abbreviations <- function(pathway_names) {
  
  # Define common abbreviations
  pathway_abbr <- c(
    "Inflammatory response" = "Inflammation",
    "Interferon gamma response" = "IFN-γ response",
    "Fatty acid metabolism" = "FA metabolism",
    "Oxidative phosphorylation" = "OXPHOS",
    "Epithelial mesenchymal transition" = "EMT",
    "Tnfa signaling via nfkb" = "TNFα via NF-κB",
    "Il6 jak stat3 signaling" = "IL6-JAK-STAT3",
    "Tgf beta signaling" = "TGF-β",
    "Wnt beta catenin signaling" = "WNT/β-catenin",
    "G2m checkpoint" = "G2/M checkpoint",
    "Myc targets v1" = "MYC targets"
  )
  
  # Add any missing pathways with their original names
  missing_pathways <- setdiff(pathway_names, names(pathway_abbr))
  if (length(missing_pathways) > 0) {
    missing_abbr <- setNames(missing_pathways, missing_pathways)
    pathway_abbr <- c(pathway_abbr, missing_abbr)
  }
  
  return(pathway_abbr)
}

#' Filter Enrichment Results
#' 
#' Filter enrichment results based on various criteria
#' 
#' @param enrichment_data Data frame with enrichment results
#' @param p_cutoff P-value cutoff
#' @param min_gene_count Minimum number of genes in pathway
#' @param pathway_types Vector of pathway types to include
#' @return Filtered data frame
#' @export
filter_enrichment_results <- function(enrichment_data,
                                     p_cutoff = 0.05,
                                     min_gene_count = 5,
                                     pathway_types = NULL) {
  
  filtered_data <- enrichment_data %>%
    filter(p.adjust <= p_cutoff,
           Count >= min_gene_count)
  
  if (!is.null(pathway_types)) {
    filtered_data <- filtered_data %>%
      filter(type %in% pathway_types)
  }
  
  return(filtered_data)
}

# ==============================================================================
# PATHWAY ANALYSIS HELPERS
# ==============================================================================

#' Compare Pathways Across Conditions
#' 
#' Compare pathway enrichment across different conditions
#' 
#' @param enrichment_results List of enrichment results
#' @param pathway_name Specific pathway to compare
#' @param clusters Vector of clusters to include
#' @return Data frame with comparison results
#' @export
compare_pathways_across_conditions <- function(enrichment_results,
                                              pathway_name = NULL,
                                              clusters = NULL) {
  
  comparison_data <- list()
  
  for (condition in names(enrichment_results)) {
    for (cluster in names(enrichment_results[[condition]])) {
      if (!is.null(enrichment_results[[condition]][[cluster]])) {
        
        cluster_data <- enrichment_results[[condition]][[cluster]]$significant_paths %>%
          mutate(condition = condition,
                 cluster = cluster)
        
        if (!is.null(pathway_name)) {
          cluster_data <- cluster_data %>%
            filter(pathway == pathway_name)
        }
        
        if (!is.null(clusters)) {
          cluster_data <- cluster_data %>%
            filter(cluster %in% clusters)
        }
        
        comparison_data[[paste(condition, cluster, sep = "_")]] <- cluster_data
      }
    }
  }
  
  combined_data <- do.call(rbind, comparison_data)
  rownames(combined_data) <- NULL
  
  return(combined_data)
}

#' Get Top Pathways
#' 
#' Extract top pathways from enrichment results
#' 
#' @param enrichment_results List of enrichment results
#' @param n_top Number of top pathways to return
#' @param metric Metric to rank by ("p.adjust", "Count", "GeneRatio")
#' @return Vector of top pathway names
#' @export
get_top_pathways <- function(enrichment_results,
                            n_top = 20,
                            metric = "p.adjust") {
  
  all_pathways <- list()
  
  for (condition in names(enrichment_results)) {
    for (cluster in names(enrichment_results[[condition]])) {
      if (!is.null(enrichment_results[[condition]][[cluster]])) {
        
        cluster_pathways <- enrichment_results[[condition]][[cluster]]$significant_paths %>%
          mutate(condition = condition, cluster = cluster)
        
        all_pathways[[paste(condition, cluster, sep = "_")]] <- cluster_pathways
      }
    }
  }
  
  combined_pathways <- do.call(rbind, all_pathways)
  
  if (metric == "p.adjust") {
    top_pathways <- combined_pathways %>%
      group_by(pathway) %>%
      summarise(best_pval = min(p.adjust), .groups = "drop") %>%
      arrange(best_pval) %>%
      slice_head(n = n_top) %>%
      pull(pathway)
  } else if (metric == "Count") {
    top_pathways <- combined_pathways %>%
      group_by(pathway) %>%
      summarise(max_count = max(Count), .groups = "drop") %>%
      arrange(desc(max_count)) %>%
      slice_head(n = n_top) %>%
      pull(pathway)
  } else if (metric == "GeneRatio") {
    top_pathways <- combined_pathways %>%
      group_by(pathway) %>%
      summarise(max_ratio = max(GeneRatio), .groups = "drop") %>%
      arrange(desc(max_ratio)) %>%
      slice_head(n = n_top) %>%
      pull(pathway)
  }
  
  return(top_pathways)
}

#' Create Pathway Summary Report
#' 
#' Generate a comprehensive summary report of pathway analysis
#' 
#' @param enrichment_results List of enrichment results
#' @param output_file Output file path for report
#' @return Summary statistics
#' @export
create_pathway_summary_report <- function(enrichment_results,
                                         output_file = "pathway_summary.txt") {
  
  summary_stats <- list()
  
  # Overall statistics
  total_conditions <- length(enrichment_results)
  total_clusters <- sum(sapply(enrichment_results, length))
  
  # Pathway statistics
  all_pathways <- compare_pathways_across_conditions(enrichment_results)
  
  if (nrow(all_pathways) > 0) {
    unique_pathways <- length(unique(all_pathways$pathway))
    pathway_types <- table(all_pathways$type)
    
    # Top pathways by frequency
    pathway_frequency <- all_pathways %>%
      count(pathway, sort = TRUE) %>%
      slice_head(n = 10)
    
    # Summary by condition
    condition_summary <- all_pathways %>%
      group_by(condition) %>%
      summarise(
        n_pathways = n(),
        n_unique_pathways = n_distinct(pathway),
        avg_pval = mean(p.adjust),
        .groups = "drop"
      )
    
    # Summary by cluster
    cluster_summary <- all_pathways %>%
      group_by(cluster) %>%
      summarise(
        n_pathways = n(),
        n_unique_pathways = n_distinct(pathway),
        avg_pval = mean(p.adjust),
        .groups = "drop"
      )
    
    summary_stats <- list(
      total_conditions = total_conditions,
      total_clusters = total_clusters,
      unique_pathways = unique_pathways,
      pathway_types = pathway_types,
      top_pathways = pathway_frequency,
      condition_summary = condition_summary,
      cluster_summary = cluster_summary
    )
    
    # Write report
    sink(output_file)
    cat("=== PATHWAY ENRICHMENT ANALYSIS SUMMARY ===\n\n")
    cat(sprintf("Total conditions analyzed: %d\n", total_conditions))
    cat(sprintf("Total clusters analyzed: %d\n", total_clusters))
    cat(sprintf("Unique significant pathways found: %d\n\n", unique_pathways))
    
    cat("Pathway types distribution:\n")
    print(pathway_types)
    cat("\n")
    
    cat("Top 10 most frequent pathways:\n")
    print(pathway_frequency)
    cat("\n")
    
    cat("Summary by condition:\n")
    print(condition_summary)
    cat("\n")
    
    cat("Summary by cluster:\n")
    print(cluster_summary)
    cat("\n")
    
    sink()
    
    cat(sprintf("Summary report written to: %s\n", output_file))
  } else {
    cat("No significant pathways found in the analysis.\n")
  }
  
  return(summary_stats)
}