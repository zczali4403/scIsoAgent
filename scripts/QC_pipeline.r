# ------------------------------------------------------------------------------
# Script: QC_pipeline.r
# Purpose: Staged QC pipeline for pre-filtering, PCA, clustering, and doublet removal.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(DoubletFinder)
})

# Parse command-line arguments.
# Brief: Parse command-line arguments written as `--key=value` pairs.
# Args:
#   args: Raw trailing command-line arguments passed to the script.
# Returns: Named list/dictionary of parsed command-line options.
parse_args <- function(args) {
  parsed <- list()

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop(
        paste(
          "Invalid argument:", arg,
          "\nExpected format: --key=value"
        )
      )
    }

    key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- key_value[1]
    value <- paste(key_value[-1], collapse = "=")
    parsed[[key]] <- value
  }

  parsed
}

# Brief: Check that required command-line options were supplied before running a stage.
# Args:
#   parsed: Named argument list produced by the script argument parser.
#   required_keys: Argument names that must be present before continuing.
# Returns: No return value; raises an error when any required option is missing.
require_args <- function(parsed, required_keys) {
  missing_keys <- required_keys[!required_keys %in% names(parsed)]
  if (length(missing_keys) > 0) {
    stop(
      paste(
        "Missing required arguments:",
        paste(missing_keys, collapse = ", ")
      )
    )
  }
}

# Brief: Create an output directory when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the target directory exists; no analysis object is returned.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Brief: Add mitochondrial-percentage metadata when it is absent from a Seurat object.
# Args:
#   seurat_object: Seurat object being processed or plotted.
# Returns: Seurat object with percent.mt metadata present or newly calculated.
ensure_percent_mt <- function(seurat_object) {
  if (!"percent.mt" %in% colnames(seurat_object@meta.data)) {
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  }
  seurat_object
}

# Brief: Write a ggplot-compatible plot object to a PDF file.
# Args:
#   plot_obj: Plot object to save.
#   file_path: Filesystem path for an input or output file.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Writes the provided plot object to a PDF file.
save_plot_pdf <- function(plot_obj, file_path, width = 7, height = 6) {
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Render a data frame as a PDF table.
# Args:
#   df: Data frame to render as a table in a PDF report.
#   file_path: Filesystem path for an input or output file.
#   title: Plot title text.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Writes the provided data frame as a PDF table.
save_table_pdf <- function(df, file_path, title = NULL, width = 8, height = 5) {
  tbl <- tableGrob(df, rows = NULL)
  pdf(file = file_path, width = width, height = height, bg = "white")
  grid.newpage()
  if (!is.null(title)) {
    grid.text(title, y = unit(0.97, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
    grid.draw(editGrob(tbl, vp = viewport(y = 0.45, height = 0.82)))
  } else {
    grid.draw(tbl)
  }
  dev.off()
}

# Brief: Create a pre-filter scatter plot of total counts versus detected features.
# Args:
#   seurat_object: Seurat object being processed or plotted.
# Returns: ggplot/plot object for scatter plot.
build_scatter_plot <- function(seurat_object) {
  md <- seurat_object@meta.data
  ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(alpha = 0.35, size = 0.4, color = "#3B82F6") +
    geom_smooth(method = "lm", se = FALSE, color = "#D9485F", linewidth = 0.7) +
    labs(
      title = "Reads vs Unique Genes per Cell BEFORE Filtering",
      x = "nCount_RNA",
      y = "nFeature_RNA"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

# Brief: Create violin plots for QC metrics or isoform-level count summaries.
# Args:
#   seurat_object: Seurat object being processed or plotted.
#   title_text: Plot title text.
# Returns: ggplot/plot object for violin plot.
build_violin_plot <- function(seurat_object, title_text) {
  md <- seurat_object@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt"), drop = FALSE]
  long_df <- stack(md)
  colnames(long_df) <- c("value", "metric")

  ggplot(long_df, aes(x = metric, y = value, fill = metric)) +
    geom_violin(scale = "width", trim = TRUE, color = "#374151", linewidth = 0.2) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    labs(title = title_text, x = NULL, y = "Value") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# Brief: Create a PCA elbow plot for selecting dimensionality.
# Args:
#   seurat_object: Seurat object being processed or plotted.
# Returns: ggplot/plot object for elbow plot.
build_elbow_plot <- function(seurat_object) {
  ElbowPlot(seurat_object, ndims = 50) +
    ggtitle("Elbow Plot") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

# Brief: Create an embedding plot colored by a metadata column.
# Args:
#   seurat_object: Seurat object being processed or plotted.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   title_text: Plot title text.
# Returns: ggplot/plot object for cell dim plot.
build_cell_dim_plot <- function(seurat_object, group_by, title_text) {
  CellDimPlot(
    srt = seurat_object,
    group.by = group_by,
    reduction = "UMAP",
    theme_use = "theme_blank",
    title = title_text
  )
}

# Brief: Create an embedding plot colored by one feature or QC metric.
# Args:
#   seurat_object: Seurat object being processed or plotted.
#   feature_name: Single feature or QC metric to display on an embedding.
#   title_text: Plot title text.
# Returns: ggplot/plot object for feature dim plot.
build_feature_dim_plot <- function(seurat_object, feature_name, title_text) {
  FeatureDimPlot(
    srt = seurat_object,
    features = feature_name,
    reduction = "UMAP",
    theme_use = "theme_blank"
  ) +
    ggtitle(title_text)
}

# Brief: Apply QC thresholds and run normalization, variable-feature selection, scaling, and PCA.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   min_features: Minimum detected features required to keep a cell or create an object.
#   max_features: Maximum detected features allowed for QC filtering.
#   min_counts: Minimum total counts required for QC filtering.
#   max_counts: Maximum total counts allowed for QC filtering.
#   mt_cutoff: Maximum mitochondrial percentage allowed for QC filtering.
# Returns: Filtered Seurat object with normalized data, variable features, scaling, and PCA.
prepare_filtered_pca_object <- function(input_rds,
                                        min_features,
                                        max_features,
                                        min_counts,
                                        max_counts,
                                        mt_cutoff) {
  seurat_object <- readRDS(input_rds)
  seurat_object <- ensure_percent_mt(seurat_object)

  filt_seurat_object <- subset(
    seurat_object,
    subset =
      nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < mt_cutoff &
      nCount_RNA < max_counts &
      nCount_RNA > min_counts
  )

  filt_seurat_object <- NormalizeData(filt_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  filt_seurat_object <- FindVariableFeatures(filt_seurat_object, selection.method = "vst", nfeatures = 2000)
  filt_seurat_object <- ScaleData(filt_seurat_object, features = rownames(filt_seurat_object))
  filt_seurat_object <- RunPCA(filt_seurat_object, features = VariableFeatures(object = filt_seurat_object))

  filt_seurat_object
}

# Brief: Build a filtered PCA object and add neighbors, clusters, and UMAP coordinates.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   min_features: Minimum detected features required to keep a cell or create an object.
#   max_features: Maximum detected features allowed for QC filtering.
#   min_counts: Minimum total counts required for QC filtering.
#   max_counts: Maximum total counts allowed for QC filtering.
#   mt_cutoff: Maximum mitochondrial percentage allowed for QC filtering.
#   npc: Number of PCA dimensions to use.
#   cluster_res: Clustering resolution used by Seurat FindClusters.
# Returns: Filtered Seurat object with neighbors, clusters, and UMAP coordinates.
prepare_clustered_object <- function(input_rds,
                                     min_features,
                                     max_features,
                                     min_counts,
                                     max_counts,
                                     mt_cutoff,
                                     npc,
                                     cluster_res) {
  filt_seurat_object <- prepare_filtered_pca_object(
    input_rds = input_rds,
    min_features = min_features,
    max_features = max_features,
    min_counts = min_counts,
    max_counts = max_counts,
    mt_cutoff = mt_cutoff
  )

  filt_seurat_object <- FindNeighbors(filt_seurat_object, dims = 1:npc)
  filt_seurat_object <- FindClusters(filt_seurat_object, resolution = cluster_res)
  filt_seurat_object <- RunUMAP(filt_seurat_object, dims = 1:npc)

  filt_seurat_object
}

# Brief: Generate QC plots before applying filtering thresholds.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   output_dir: Directory where output files should be written.
#   sample_id: Sample identifier used in metadata and output filenames.
# Returns: Writes pre-filter QC scatter and violin PDFs.
pre_filter_stage <- function(input_rds, output_dir, sample_id) {
  seurat_object <- readRDS(input_rds)
  seurat_object <- ensure_percent_mt(seurat_object)
  ensure_dir(output_dir)

  scatter_plot <- build_scatter_plot(seurat_object)
  violin_plot <- build_violin_plot(seurat_object, "QC plots (gene level) BEFORE Filtering")

  save_plot_pdf(scatter_plot, file.path(output_dir, "nCount_nFeature_corr.pdf"))
  save_plot_pdf(violin_plot, file.path(output_dir, "QC_metrics_violin.pdf"), width = 12, height = 5)
}

# Brief: Apply user-confirmed QC thresholds and generate post-filter/PCA diagnostics.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   output_dir: Directory where output files should be written.
#   sample_id: Sample identifier used in metadata and output filenames.
#   min_features: Minimum detected features required to keep a cell or create an object.
#   max_features: Maximum detected features allowed for QC filtering.
#   min_counts: Minimum total counts required for QC filtering.
#   max_counts: Maximum total counts allowed for QC filtering.
#   mt_cutoff: Maximum mitochondrial percentage allowed for QC filtering.
# Returns: Writes post-filter violin plot, PCA elbow plot, and filter parameter CSV.
filter_pca_stage <- function(input_rds,
                             output_dir,
                             sample_id,
                             min_features,
                             max_features,
                             min_counts,
                             max_counts,
                             mt_cutoff) {
  ensure_dir(output_dir)

  seurat_object <- readRDS(input_rds)
  seurat_object <- ensure_percent_mt(seurat_object)
  filt_seurat_object <- prepare_filtered_pca_object(
    input_rds = input_rds,
    min_features = min_features,
    max_features = max_features,
    min_counts = min_counts,
    max_counts = max_counts,
    mt_cutoff = mt_cutoff
  )

  violin_after <- build_violin_plot(filt_seurat_object, "QC metrics gene level AFTER Filtering")
  save_plot_pdf(violin_after, file.path(output_dir, "QC_metrics_after_qc_violin.pdf"), width = 12, height = 5)

  elbow_plot <- build_elbow_plot(filt_seurat_object)
  save_plot_pdf(elbow_plot, file.path(output_dir, "Elbowplot.pdf"))

  filter_params <- data.frame(
    parameter = c("min_features", "max_features", "min_counts", "max_counts", "mt_cutoff"),
    value = c(min_features, max_features, min_counts, max_counts, mt_cutoff)
  )
  write.csv(filter_params, file.path(output_dir, "filter_params.csv"), row.names = FALSE)
}

# Brief: Run clustering and UMAP using user-confirmed PCA dimensions and resolution.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   output_dir: Directory where output files should be written.
#   sample_id: Sample identifier used in metadata and output filenames.
#   min_features: Minimum detected features required to keep a cell or create an object.
#   max_features: Maximum detected features allowed for QC filtering.
#   min_counts: Minimum total counts required for QC filtering.
#   max_counts: Maximum total counts allowed for QC filtering.
#   mt_cutoff: Maximum mitochondrial percentage allowed for QC filtering.
#   npc: Number of PCA dimensions to use.
#   cluster_res: Clustering resolution used by Seurat FindClusters.
# Returns: Runs clustering/UMAP for validation; no file is written by this stage function itself.
cluster_umap_stage <- function(input_rds,
                               output_dir,
                               sample_id,
                               min_features,
                               max_features,
                               min_counts,
                               max_counts,
                               mt_cutoff,
                               npc,
                               cluster_res) {
  ensure_dir(output_dir)

  filt_seurat_object <- prepare_clustered_object(
    input_rds = input_rds,
    min_features = min_features,
    max_features = max_features,
    min_counts = min_counts,
    max_counts = max_counts,
    mt_cutoff = mt_cutoff,
    npc = npc,
    cluster_res = cluster_res
  )

}

# Brief: Run DoubletFinder, keep singlets, and save final QC outputs.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   output_dir: Directory where output files should be written.
#   sample_id: Sample identifier used in metadata and output filenames.
#   min_features: Minimum detected features required to keep a cell or create an object.
#   max_features: Maximum detected features allowed for QC filtering.
#   min_counts: Minimum total counts required for QC filtering.
#   max_counts: Maximum total counts allowed for QC filtering.
#   mt_cutoff: Maximum mitochondrial percentage allowed for QC filtering.
#   npc: Number of PCA dimensions to use.
#   cluster_res: Clustering resolution used by Seurat FindClusters.
#   doublet_rate: Expected doublet rate used by DoubletFinder.
# Returns: Writes doublet diagnostics, QC summaries, and final singlet/doublet Seurat RDS files.
doublet_stage <- function(input_rds,
                          output_dir,
                          sample_id,
                          min_features,
                          max_features,
                          min_counts,
                          max_counts,
                          mt_cutoff,
                          npc,
                          cluster_res,
                          doublet_rate = 0.05) {
  ensure_dir(output_dir)

  filt_seurat_object <- prepare_clustered_object(
    input_rds = input_rds,
    min_features = min_features,
    max_features = max_features,
    min_counts = min_counts,
    max_counts = max_counts,
    mt_cutoff = mt_cutoff,
    npc = npc,
    cluster_res = cluster_res
  )

  sweep.res.list_pbmc <- paramSweep(filt_seurat_object, PCs = 1:npc, sct = FALSE)
  sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
  bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
  pK <- bcmvn_pbmc %>%
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))

  annotations <- filt_seurat_object@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate * nrow(filt_seurat_object@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  filt_seurat_object <- doubletFinder(
    filt_seurat_object,
    PCs = 1:npc,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  colnames(filt_seurat_object@meta.data) <- sub(
    "DF.classifications_.*$",
    "DF.classifications",
    colnames(filt_seurat_object@meta.data)
  )

  statsDoublets <- filt_seurat_object@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(
      Median_nCount_RNA = median(nCount_RNA),
      Median_nFeature_RNA = median(nFeature_RNA),
      Count = n(),
      .groups = "drop"
    )

  filt_seurat_object_doublets <- filt_seurat_object
  filt_seurat_object <- subset(filt_seurat_object, subset = DF.classifications == "Singlet")

  p_clusters <- build_cell_dim_plot(filt_seurat_object, "seurat_clusters", "Seurat Clusters UMAP")
  p_ncount <- build_feature_dim_plot(filt_seurat_object, "nCount_RNA", "nCount_RNA UMAP")
  p_nfeature <- build_feature_dim_plot(filt_seurat_object, "nFeature_RNA", "nFeature_RNA UMAP")
  p_doublet <- build_cell_dim_plot(filt_seurat_object_doublets, "DF.classifications", "Doublet Classification UMAP")

  save_plot_pdf(p_clusters, file.path(output_dir, "seurat_clusters_umap.pdf"))
  save_plot_pdf(p_ncount, file.path(output_dir, "nCount_RNA_umap.pdf"))
  save_plot_pdf(p_nfeature, file.path(output_dir, "nFeature_RNA_umap.pdf"))
  save_plot_pdf(p_doublet, file.path(output_dir, "DF_classification_umap.pdf"))

  save_table_pdf(statsDoublets, file.path(output_dir, "stats_Doublets.pdf"), "Doublet Summary")

  stats_summary <- data.frame(
    metric = c(
      "Sample ID",
      "Cells_before_filter",
      "Cells_after_filter",
      "Median Feature per Cell",
      "Median Reads per Cell",
      "NPCs",
      "Doublet Rate",
      "Chosen pK"
    ),
    value = c(
      sample_id,
      ncol(filt_seurat_object_doublets),
      ncol(filt_seurat_object),
      median(filt_seurat_object$nFeature_RNA),
      median(filt_seurat_object$nCount_RNA),
      npc,
      doublet_rate,
      pK
    )
  )
  save_table_pdf(stats_summary, file.path(output_dir, "stats_summary.pdf"), "QC Summary")

  saveRDS(filt_seurat_object_doublets, file.path(output_dir, paste0(sample_id, "_with_doublets.rds")))
  saveRDS(filt_seurat_object, file.path(output_dir, paste0(sample_id, "_singlets.rds")))

  write.csv(bcmvn_pbmc, file.path(output_dir, "doubletFinder_pK_scan.csv"), row.names = FALSE)
}

# Main workflow entrypoint.
# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_args(args, c("stage", "input_rds", "output_dir", "sample_id"))

  stage <- args$stage
  input_rds <- args$input_rds
  output_dir <- args$output_dir
  sample_id <- args$sample_id

  if (stage == "pre-filter") {
    pre_filter_stage(input_rds, output_dir, sample_id)
    return(invisible(NULL))
  }

  if (stage == "filter-pca") {
    require_args(args, c("min_features", "max_features", "min_counts", "max_counts", "mt_cutoff"))
    filter_pca_stage(
      input_rds = input_rds,
      output_dir = output_dir,
      sample_id = sample_id,
      min_features = as.numeric(args$min_features),
      max_features = as.numeric(args$max_features),
      min_counts = as.numeric(args$min_counts),
      max_counts = as.numeric(args$max_counts),
      mt_cutoff = as.numeric(args$mt_cutoff)
    )
    return(invisible(NULL))
  }

  if (stage == "cluster-umap") {
    require_args(args, c("min_features", "max_features", "min_counts", "max_counts", "mt_cutoff", "npc", "cluster_res"))
    cluster_umap_stage(
      input_rds = input_rds,
      output_dir = output_dir,
      sample_id = sample_id,
      min_features = as.numeric(args$min_features),
      max_features = as.numeric(args$max_features),
      min_counts = as.numeric(args$min_counts),
      max_counts = as.numeric(args$max_counts),
      mt_cutoff = as.numeric(args$mt_cutoff),
      npc = as.integer(args$npc),
      cluster_res = as.numeric(args$cluster_res)
    )
    return(invisible(NULL))
  }

  if (stage == "doublet") {
    require_args(args, c("min_features", "max_features", "min_counts", "max_counts", "mt_cutoff", "npc", "cluster_res"))
    doublet_rate <- if ("doublet_rate" %in% names(args)) as.numeric(args$doublet_rate) else 0.05
    doublet_stage(
      input_rds = input_rds,
      output_dir = output_dir,
      sample_id = sample_id,
      min_features = as.numeric(args$min_features),
      max_features = as.numeric(args$max_features),
      min_counts = as.numeric(args$min_counts),
      max_counts = as.numeric(args$max_counts),
      mt_cutoff = as.numeric(args$mt_cutoff),
      npc = as.integer(args$npc),
      cluster_res = as.numeric(args$cluster_res),
      doublet_rate = doublet_rate
    )
    return(invisible(NULL))
  }

  stop("Unsupported --stage. Use one of: pre-filter, filter-pca, cluster-umap, doublet")
}

main()
