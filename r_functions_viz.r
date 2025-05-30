#' ScellPortal: Visualization Functions
#' 
#' Comprehensive plotting functions for single-cell analysis
#' 
#' @author Your Name

# ==============================================================================
# CORE PLOTTING FUNCTIONS
# ==============================================================================

#' Create UMAP Plot
#' 
#' Generate customizable UMAP visualization
#' 
#' @param obj Seurat object
#' @param group_by Metadata column to color by
#' @param reduction Reduction to use for coordinates
#' @param colors Custom color palette
#' @param pt_size Point size
#' @param add_labels Whether to add cluster labels
#' @param title Plot title
#' @return ggplot object
#' @export
create_umap_plot <- function(obj, 
                            group_by = "seurat_clusters",
                            reduction = "umap",
                            colors = NULL,
                            pt_size = 0.5,
                            add_labels = FALSE,
                            title = NULL) {
  
  p <- DimPlot(obj, 
               reduction = reduction, 
               group.by = group_by,
               cols = colors,
               pt.size = pt_size) +
    labs(title = title) +
    scellportal_theme()
  
  if (add_labels) {
    # Add cluster labels at centroids
    df <- data.frame(
      UMAP_1 = obj@reductions[[reduction]]@cell.embeddings[, 1],
      UMAP_2 = obj@reductions[[reduction]]@cell.embeddings[, 2],
      cluster = obj@meta.data[[group_by]]
    ) %>%
      group_by(cluster) %>%
      summarise(
        center_x = median(UMAP_1),
        center_y = median(UMAP_2),
        .groups = 'drop'
      )
    
    p <- p + 
      geom_point(data = df, aes(center_x, center_y), 
                 size = 10, color = "grey90", alpha = 0.6, inherit.aes = FALSE) +
      geom_text(data = df, aes(center_x, center_y, label = cluster), 
                color = "black", size = 5, fontface = "bold", inherit.aes = FALSE)
  }
  
  return(p)
}

#' Create Violin Plots
#' 
#' Generate violin plots for gene expression
#' 
#' @param obj Seurat object
#' @param genes Vector of genes to plot
#' @param group_by Grouping variable
#' @param ncol Number of columns in facet
#' @param colors Color palette
#' @param add_stats Whether to add statistical comparisons
#' @return ggplot object
#' @export
create_violin_plots <- function(obj, 
                               genes, 
                               group_by = "seurat_clusters",
                               ncol = 3,
                               colors = NULL,
                               add_stats = FALSE) {
  
  # Convert gene symbols to IDs if needed
  if (exists("bm5")) {
    gene_ids <- bm5$ensembl_gene_id[match(genes, bm5$hgnc_symbol)]
    gene_ids <- gene_ids[!is.na(gene_ids)]
  } else {
    gene_ids <- genes
  }
  
  if (add_stats && requireNamespace("ggpubr", quietly = TRUE)) {
    # Create violin plot with statistics
    plot_data <- FetchData(obj, vars = c(group_by, gene_ids))
    colnames(plot_data)[1] <- "Group"
    colnames(plot_data)[2:ncol(plot_data)] <- genes[1:length(gene_ids)]
    
    # Reshape data
    plot_data_long <- plot_data %>%
      tidyr::pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")
    
    # Get comparison groups
    groups <- levels(factor(plot_data_long$Group))
    comparisons <- combn(groups, 2, simplify = FALSE)
    
    p <- ggpubr::ggviolin(plot_data_long, x = "Group", y = "Expression", 
                          fill = "Group", palette = colors) +
      ggpubr::stat_compare_means(comparisons = comparisons, 
                                label = "p.signif", vjust = -0.5) +
      facet_wrap(~Gene, ncol = ncol, scales = "free_y") +
      scellportal_theme() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  } else {
    # Standard Seurat violin plot
    p <- VlnPlot(obj, 
                 features = gene_ids, 
                 group.by = group_by,
                 cols = colors,
                 ncol = ncol,
                 pt.size = 0) +
      scellportal_theme()
    
    # Update labels with gene symbols
    if (exists("bm5")) {
      p$data$features.plot <- factor(
        bm5$hgnc_symbol[match(p$data$features.plot, bm5$ensembl_gene_id)],
        levels = genes
      )
    }
  }
  
  return(p)
}

#' Create Density Plots
#' 
#' Generate density plots using Nebulosa
#' 
#' @param obj Seurat object
#' @param genes Vector of genes to plot
#' @param reduction Reduction to use
#' @param pt_size Point size
#' @param palette Color palette
#' @return List of ggplot objects
#' @export
create_density_plots <- function(obj, 
                                genes,
                                reduction = "umap",
                                pt_size = 0.4,
                                palette = "magma") {
  
  if (!requireNamespace("Nebulosa", quietly = TRUE)) {
    stop("Nebulosa package required for density plots")
  }
  
  # Convert gene symbols to IDs if needed
  if (exists("bm5")) {
    gene_ids <- bm5$ensembl_gene_id[match(genes, bm5$hgnc_symbol)]
  } else {
    gene_ids <- genes
  }
  
  plots <- list()
  
  for (i in seq_along(genes)) {
    p <- Nebulosa::plot_density(
      obj,
      features = gene_ids[i],
      slot = "data",
      reduction = reduction,
      size = pt_size,
      pal = palette
    ) +
      ggtitle(genes[i]) +
      density_plot_theme()
    
    plots[[genes[i]]] <- p
  }
  
  return(plots)
}

#' Create Pie Charts
#' 
#' Generate pie charts for cell composition analysis
#' 
#' @param obj Seurat object
#' @param group_by Primary grouping variable
#' @param split_by Secondary grouping variable
#' @param colors Color palette
#' @return ggplot object
#' @export
create_pie_charts <- function(obj, 
                             group_by = "cell_type",
                             split_by = "condition",
                             colors = NULL) {
  
  # Calculate cell counts
  cell_counts <- table(obj@meta.data[[group_by]], obj@meta.data[[split_by]])
  data_long <- as.data.frame(cell_counts)
  colnames(data_long) <- c("cell_type", "condition", "Count")
  
  # Create pie chart function
  create_single_pie <- function(condition_name) {
    condition_data <- data_long %>%
      filter(condition == condition_name) %>%
      mutate(
        Percentage = Count / sum(Count) * 100,
        ymax = cumsum(Percentage),
        ymin = c(0, head(cumsum(Percentage), n = -1)),
        labelPosition = (ymax + ymin) / 2,
        label = paste0(round(Percentage, 1), "%")
      )
    
    ggplot(condition_data, aes(ymax = ymax, ymin = ymin, 
                              xmax = 4, xmin = 2, fill = cell_type)) +
      geom_rect(color = "white", size = 0.5) +
      coord_polar(theta = "y") +
      xlim(c(0, 4)) +
      theme_void() +
      scale_fill_manual(values = colors) +
      geom_text(aes(label = label, x = 3, y = labelPosition),
                size = 3.5, color = "black", fontface = "bold") +
      labs(title = condition_name, fill = "Cell Types") +
      pie_chart_theme()
  }
  
  # Create pie charts for each condition
  conditions <- unique(data_long$condition)
  pie_plots <- map(conditions, create_single_pie)
  names(pie_plots) <- conditions
  
  return(pie_plots)
}

#' Create Heatmap
#' 
#' Generate expression heatmap with custom annotations
#' 
#' @param obj Seurat object
#' @param genes Vector of genes to plot
#' @param group_by Grouping variable
#' @param gene_groups Named list of gene groups
#' @param colors Color scheme
#' @return ggplot object
#' @export
create_heatmap <- function(obj, 
                          genes,
                          group_by = "seurat_clusters",
                          gene_groups = NULL,
                          colors = c("blue", "white", "red")) {
  
  # Convert genes to IDs if needed
  if (exists("bm5")) {
    gene_ids <- bm5$ensembl_gene_id[match(genes, bm5$hgnc_symbol)]
  } else {
    gene_ids <- genes
  }
  
  # Create dot plot data
  dot_plot <- DotPlot(obj, features = gene_ids, group.by = group_by)
  plot_data <- dot_plot$data
  
  # Add gene groups if provided
  if (!is.null(gene_groups)) {
    plot_data$gene_group <- NA
    for (group_name in names(gene_groups)) {
      group_genes <- gene_groups[[group_name]]
      if (exists("bm5")) {
        group_gene_ids <- bm5$ensembl_gene_id[match(group_genes, bm5$hgnc_symbol)]
      } else {
        group_gene_ids <- group_genes
      }
      plot_data$gene_group[plot_data$features.plot %in% group_gene_ids] <- group_name
    }
  }
  
  # Convert back to gene symbols for plotting
  if (exists("bm5")) {
    plot_data$features.plot <- bm5$hgnc_symbol[match(plot_data$features.plot, bm5$ensembl_gene_id)]
  }
  
  # Create heatmap
  p <- ggplot(plot_data, aes(x = id, y = features.plot, fill = avg.exp.scaled)) +
    geom_tile() +
    scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3],
                        name = "Scaled\nExpression") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black")
    )
  
  # Add gene group separators if provided
  if (!is.null(gene_groups)) {
    # Add horizontal lines between gene groups
    group_breaks <- cumsum(sapply(gene_groups, length)) + 0.5
    p <- p + geom_hline(yintercept = group_breaks[-length(group_breaks)], 
                       color = "black", size = 1)
  }
  
  return(p)
}

# ==============================================================================
# DIFFERENTIAL EXPRESSION VISUALIZATION
# ==============================================================================

#' Create Volcano Plot
#' 
#' Generate volcano plot for differential expression results
#' 
#' @param deg_results Data frame with DE results
#' @param genes_to_label Vector of genes to label
#' @param fc_cutoff Log2 fold change cutoff
#' @param p_cutoff Adjusted p-value cutoff
#' @param colors Named vector of colors for different categories
#' @return ggplot object
#' @export
create_volcano_plot <- function(deg_results,
                               genes_to_label = NULL,
                               fc_cutoff = 1,
                               p_cutoff = 0.05,
                               colors = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "gray40")) {
  
  # Add significance categories
  deg_results <- deg_results %>%
    mutate(
      significance = case_when(
        avg_log2FC > fc_cutoff & p_val_adj < p_cutoff ~ "UP",
        avg_log2FC < -fc_cutoff & p_val_adj < p_cutoff ~ "DOWN",
        TRUE ~ "NS"
      ),
      neg_log10_p = -log10(p_val_adj)
    )
  
  # Create base plot
  p <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_p, color = significance)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = colors, 
                      breaks = c("DOWN", "UP"),
                      labels = c("Down-regulated", "Up-regulated")) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), 
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(p_cutoff), 
               linetype = "dashed", color = "gray50") +
    labs(x = expression(paste("Log"[2], " fold change")),
         y = expression("-Log"[10]*" p-value"),
         color = NULL) +
    theme_minimal() +
    theme(
      legend.position = "top",
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black")
    )
  
  # Add gene labels if provided
  if (!is.null(genes_to_label)) {
    label_data <- deg_results %>%
      filter(gene %in% genes_to_label)
    
    if (nrow(label_data) > 0) {
      p <- p + 
        ggrepel::geom_text_repel(
          data = label_data,
          aes(label = gene),
          color = "black",
          fontface = "bold",
          size = 3.5,
          max.overlaps = 20,
          box.padding = 0.5
        )
    }
  }
  
  return(p)
}

#' Create Waterfall Plot
#' 
#' Generate waterfall plot showing ranked fold changes
#' 
#' @param deg_results Data frame with DE results
#' @param n_labels Number of top/bottom genes to label
#' @return ggplot object
#' @export
create_waterfall_plot <- function(deg_results, n_labels = 10) {
  
  # Rank genes by fold change
  deg_ranked <- deg_results %>%
    arrange(avg_log2FC) %>%
    mutate(rank = row_number()) %>%
    filter(!is.na(gene))
  
  # Get top and bottom genes for labeling
  top_genes <- deg_ranked %>% 
    arrange(desc(avg_log2FC)) %>% 
    slice_head(n = n_labels)
  
  bottom_genes <- deg_ranked %>% 
    arrange(avg_log2FC) %>% 
    slice_head(n = n_labels)
  
  # Create plot
  p <- ggplot(deg_ranked, aes(x = rank, y = avg_log2FC)) +
    geom_point(aes(size = abs(avg_log2FC)), alpha = 0.7, color = "gray50") +
    geom_point(data = rbind(top_genes, bottom_genes), 
               aes(size = abs(avg_log2FC)), color = "red", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_size_continuous(range = c(0.5, 3), guide = "none") +
    labs(x = "Gene Rank", y = "Log2 Fold Change") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  # Add gene labels
  if (nrow(top_genes) > 0) {
    p <- p + 
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        nudge_x = 50,
        nudge_y = 0.2,
        color = "black",
        fontface = "bold",
        size = 3
      )
  }
  
  if (nrow(bottom_genes) > 0) {
    p <- p + 
      ggrepel::geom_text_repel(
        data = bottom_genes,
        aes(label = gene),
        nudge_x = 50,
        nudge_y = -0.2,
        color = "black",
        fontface = "bold",
        size = 3
      )
  }
  
  return(p)
}

# ==============================================================================
# PATHWAY ENRICHMENT VISUALIZATION
# ==============================================================================

#' Create Enrichment Dot Plot
#' 
#' Generate dot plot for pathway enrichment results
#' 
#' @param enrichment_data List of enrichment results
#' @param conditions Vector of condition names
#' @param top_pathways Number of top pathways to show
#' @param pathway_colors Named vector of pathway colors
#' @return ggplot object
#' @export
create_enrichment_dotplot <- function(enrichment_data,
                                     conditions = NULL,
                                     top_pathways = 20,
                                     pathway_colors = NULL) {
  
  # Combine enrichment data from all conditions
  combined_data <- map_dfr(names(enrichment_data), function(condition) {
    if (!is.null(enrichment_data[[condition]])) {
      enrichment_data[[condition]] %>%
        mutate(condition = condition)
    }
  })
  
  if (nrow(combined_data) == 0) {
    stop("No enrichment data found")
  }
  
  # Filter for top pathways and significance
  plot_data <- combined_data %>%
    filter(type == "HALLMARK", p.adjust < 0.01) %>%
    mutate(
      neg_log10_p = -log10(p.adjust),
      neg_log10_p = pmin(neg_log10_p, 20)  # Cap extreme values
    ) %>%
    group_by(pathway) %>%
    filter(max(neg_log10_p) > 2) %>%  # Only pathways with significant results
    ungroup()
  
  # Select top pathways based on maximum significance
  top_pathway_names <- plot_data %>%
    group_by(pathway) %>%
    summarise(max_sig = max(neg_log10_p), .groups = "drop") %>%
    arrange(desc(max_sig)) %>%
    slice_head(n = top_pathways) %>%
    pull(pathway)
  
  plot_data <- plot_data %>%
    filter(pathway %in% top_pathway_names) %>%
    mutate(pathway = factor(pathway, levels = rev(top_pathway_names)))
  
  # Create dot plot
  p <- ggplot(plot_data, aes(x = condition, y = pathway, 
                            fill = neg_log10_p, size = GeneRatio)) +
    geom_point(shape = 21, stroke = 0.1) +
    scale_size_continuous(range = c(2, 8), name = "Gene Ratio") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                        name = "-log10(p.adj)") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.position = "right"
    )
  
  # Add pathway colors if provided
  if (!is.null(pathway_colors)) {
    pathway_text_colors <- pathway_colors[levels(plot_data$pathway)]
    pathway_text_colors[is.na(pathway_text_colors)] <- "black"
    
    p <- p + theme(axis.text.y = element_text(color = pathway_text_colors))
  }
  
  return(p)
}

# ==============================================================================
# THEME FUNCTIONS
# ==============================================================================

#' ScellPortal Default Theme
#' 
#' Consistent theme for all plots
#' 
#' @return ggplot theme
#' @export
scellportal_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      strip.text = element_text(size = 11, face = "bold"),
      legend.key.size = unit(0.5, "cm")
    )
}

#' Density Plot Theme
#' 
#' Specialized theme for density plots
#' 
#' @return ggplot theme
#' @export
density_plot_theme <- function() {
  theme_void() +
    theme(
      text = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
}

#' Pie Chart Theme
#' 
#' Specialized theme for pie charts
#' 
#' @return ggplot theme
#' @export
pie_chart_theme <- function() {
  theme_void() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm")
    )
}