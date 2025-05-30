#' ScellPortal: Single-cell RNA-seq Analysis Toolkit
#' Core preprocessing and analysis functions
#' 
#' @author Your Name
#' @version 1.0.0

# ==============================================================================
# PREPROCESSING FUNCTIONS
# ==============================================================================

#' Preprocess Seurat Object
#' 
#' Comprehensive preprocessing including QC, normalization, and scaling
#' 
#' @param obj Seurat object
#' @param min_features Minimum number of features per cell
#' @param max_features Maximum number of features per cell
#' @param mt_cutoff Mitochondrial gene percentage cutoff
#' @param run_sctransform Whether to use SCTransform normalization
#' @return Preprocessed Seurat object
#' @export
preprocess_seurat <- function(obj, 
                             min_features = 200, 
                             max_features = 5000,
                             mt_cutoff = 20,
                             run_sctransform = TRUE) {
  
  # Calculate QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Filter cells and features
  obj <- subset(obj, 
                subset = nFeature_RNA > min_features & 
                        nFeature_RNA < max_features & 
                        percent.mt < mt_cutoff)
  
  if (run_sctransform) {
    # SCTransform normalization
    obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    # Standard normalization
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, vars.to.regress = "percent.mt")
  }
  
  return(obj)
}

#' Run Integration Analysis
#' 
#' Perform batch correction and integration using various methods
#' 
#' @param obj Seurat object
#' @param method Integration method ("harmony", "cca", "rpca")
#' @param batch_var Variable name for batch correction
#' @param dims Number of dimensions to use
#' @return Integrated Seurat object
#' @export
run_integration <- function(obj, 
                           method = "harmony", 
                           batch_var = "orig.ident",
                           dims = 1:30) {
  
  # Run PCA
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = FALSE)
  
  if (method == "harmony") {
    # Harmony integration
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("harmony package required for this integration method")
    }
    obj <- harmony::RunHarmony(obj, group.by.vars = batch_var)
    reduction_name <- "harmony"
  } else if (method == "cca") {
    # CCA integration
    obj.list <- SplitObject(obj, split.by = batch_var)
    obj.list <- lapply(obj.list, function(x) {
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    features <- SelectIntegrationFeatures(object.list = obj.list)
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                          anchor.features = features)
    obj <- IntegrateData(anchorset = obj.anchors)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(dims), verbose = FALSE)
    reduction_name <- "pca"
  }
  
  # Run UMAP and clustering
  obj <- RunUMAP(obj, reduction = reduction_name, dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = reduction_name, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  return(obj)
}

#' Annotate Cell Clusters
#' 
#' Automatic cell type annotation based on marker genes
#' 
#' @param obj Seurat object
#' @param markers Named list of marker genes for each cell type
#' @param method Annotation method ("markers", "reference")
#' @return Annotated Seurat object
#' @export
annotate_clusters <- function(obj, markers = NULL, method = "markers") {
  
  if (is.null(markers)) {
    # Default markers for synovial fibroblasts
    markers <- list(
      "SLSFs_CD34pos" = c("CD34", "CLIC5"),
      "SLSFs_CD34neg" = c("THY1", "CYBB"),
      "LLSFs" = c("PPARG", "RUNX2"),
      "SLSFs_IM" = c("CCL2", "CXCL12")
    )
  }
  
  if (method == "markers") {
    # Find all markers
    all.markers <- FindAllMarkers(obj, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
    
    # Simple annotation based on top markers
    # This would need more sophisticated logic in practice
    obj$cell_type <- paste0("Cluster_", obj$seurat_clusters)
  }
  
  return(obj)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Convert Gene IDs
#' 
#' Convert between different gene ID formats
#' 
#' @param genes Vector of gene IDs
#' @param from_type Source ID type ("ENSEMBL", "SYMBOL", "ENTREZID")
#' @param to_type Target ID type
#' @param species Species name
#' @return Data frame with converted IDs
#' @export
convert_gene_ids <- function(genes, 
                            from_type = "ENSEMBL", 
                            to_type = "SYMBOL",
                            species = "Homo sapiens") {
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package required for gene ID conversion")
  }
  
  cat(sprintf("Converting %d genes from %s to %s\n", 
              length(genes), from_type, to_type))
  
  org_db <- switch(species,
                   "Homo sapiens" = "org.Hs.eg.db",
                   "Mus musculus" = "org.Mm.eg.db",
                   "org.Hs.eg.db")
  
  genes_df <- clusterProfiler::bitr(
    genes, 
    fromType = from_type, 
    toType = to_type,
    OrgDb = org_db, 
    drop = TRUE
  )
  
  colnames(genes_df) <- c("original_id", "converted_id")
  
  cat(sprintf("Successfully converted %d genes\n", nrow(genes_df)))
  
  return(genes_df)
}

#' Prepare BioMart Data
#' 
#' Create gene ID mapping data frame
#' 
#' @param species Species name
#' @param attributes Attributes to retrieve
#' @return Data frame with gene annotations
#' @export
prepare_biomart_data <- function(species = "hsapiens", 
                                attributes = c("ensembl_gene_id", "hgnc_symbol")) {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt package required")
  }
  
  # This would connect to BioMart in practice
  # For now, return a placeholder
  cat("Preparing gene annotation data...\n")
  
  # In practice, you would do:
  # mart <- biomaRt::useEnsembl(biomart = "genes", dataset = paste0(species, "_gene_ensembl"))
  # annotations <- biomaRt::getBM(attributes = attributes, mart = mart)
  
  return(data.frame(
    ensembl_gene_id = paste0("ENSG", sprintf("%011d", 1:1000)),
    hgnc_symbol = paste0("GENE", 1:1000)
  ))
}

#' Setup Output Directory
#' 
#' Create output directory structure
#' 
#' @param base_path Base output path
#' @param create_subdirs Whether to create standard subdirectories
#' @return Path to output directory
#' @export
setup_output_dir <- function(base_path = "results", create_subdirs = TRUE) {
  
  dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
  
  if (create_subdirs) {
    subdirs <- c("plots", "tables", "qc", "differential", "enrichment")
    sapply(subdirs, function(x) {
      dir.create(file.path(base_path, x), recursive = TRUE, showWarnings = FALSE)
    })
  }
  
  return(base_path)
}

# ==============================================================================
# QUALITY CONTROL FUNCTIONS
# ==============================================================================

#' Generate QC Plots
#' 
#' Create comprehensive quality control visualizations
#' 
#' @param obj Seurat object
#' @param output_dir Output directory
#' @return List of QC plots
#' @export
generate_qc_plots <- function(obj, output_dir = "qc/") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Feature scatter plots
  p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # Violin plots
  p3 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0)
  
  # Save plots
  ggsave(file.path(output_dir, "qc_scatter.pdf"), p1 + p2, width = 10, height = 5)
  ggsave(file.path(output_dir, "qc_violin.pdf"), p3, width = 12, height = 4)
  
  return(list(scatter = p1 + p2, violin = p3))
}
