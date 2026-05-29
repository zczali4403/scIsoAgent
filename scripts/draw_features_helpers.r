# ------------------------------------------------------------------------------
# Script: draw_features_helpers.r
# Purpose: Helper utilities for iterative FeatureDimPlot visualization in persistent R sessions.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(ggplot2)
})

# Brief: Create an output directory when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the target directory exists; no analysis object is returned.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Brief: Load the Seurat object used for repeated feature plotting.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
# Returns: Loaded Seurat object for feature plotting.
read_feature_plot_object <- function(input_rds) {
  readRDS(input_rds)
}

# Brief: List assay names available in a Seurat object.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
# Returns: Character vector of assay names available in the Seurat object.
list_available_assays <- function(seurat_obj) {
  names(seurat_obj@assays)
}

# Brief: List dimensional reductions available in a Seurat object.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
# Returns: Character vector of reduction names available in the Seurat object.
list_available_reductions <- function(seurat_obj) {
  names(seurat_obj@reductions)
}

# Brief: Choose the requested assay or fall back to the active Seurat assay.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
# Returns: Resolved assay value or matching records.
resolve_assay <- function(seurat_obj,
                          assay = NULL) {
  if (is.null(assay) || assay == "") {
    return(DefaultAssay(seurat_obj))
  }

  assay_names <- list_available_assays(seurat_obj)

  if (!assay %in% assay_names) {
    stop(
      paste0(
        "Requested assay not found: ", assay,
        ". Available assays: ", paste(assay_names, collapse = ", ")
      )
    )
  }

  assay
}

# Brief: Resolve a requested reduction name with case-insensitive matching.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   reduction: Dimensional reduction name used for embedding plots.
# Returns: Resolved reduction value or matching records.
resolve_reduction <- function(seurat_obj,
                              reduction = "UMAP") {
  reduction_names <- list_available_reductions(seurat_obj)

  if (reduction %in% reduction_names) {
    return(reduction)
  }

  lower_map <- setNames(reduction_names, tolower(reduction_names))
  reduction_lower <- tolower(reduction)

  if (reduction_lower %in% names(lower_map)) {
    return(lower_map[[reduction_lower]])
  }

  stop(
    paste0(
      "Requested reduction not found: ", reduction,
      ". Available reductions: ", paste(reduction_names, collapse = ", ")
    )
  )
}

# Brief: Draw FeatureDimPlot panels for a batch of features in a persistent R session.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   reduction: Dimensional reduction name used for embedding plots.
#   theme_use: SCP plotting theme name.
# Returns: Plot object for feature batch; writes a PDF when an output path is provided.
plot_feature_batch <- function(seurat_obj,
                               features,
                               assay = NULL,
                               reduction = "UMAP",
                               theme_use = "theme_blank") {
  assay_use <- resolve_assay(seurat_obj, assay)
  reduction_use <- resolve_reduction(seurat_obj, reduction)

  DefaultAssay(seurat_obj) <- assay_use

  FeatureDimPlot(
    srt = seurat_obj,
    features = features,
    reduction = reduction_use,
    theme_use = theme_use
  )
}

# Brief: Save a feature plot produced in an interactive plotting session.
# Args:
#   plot_obj: Plot object to save.
#   output_file: Destination path for the generated table or figure.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Writes a feature plot PDF.
save_feature_plot <- function(plot_obj,
                              output_file,
                              width = 12,
                              height = 10) {
  ensure_dir(dirname(output_file))
  pdf(output_file, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Create a small table describing a feature-plot request.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   reduction: Dimensional reduction name used for embedding plots.
#   output_file: Destination path for the generated table or figure.
# Returns: ggplot/plot object for feature plot summary.
build_feature_plot_summary <- function(seurat_obj,
                                       features,
                                       assay = NULL,
                                       reduction = "UMAP",
                                       output_file = NULL) {
  data.frame(
    input_rds = NA_character_,
    assay = resolve_assay(seurat_obj, assay),
    reduction = resolve_reduction(seurat_obj, reduction),
    n_features = length(features),
    features = paste(features, collapse = ","),
    output_file = output_file %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
