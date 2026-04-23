suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(ggplot2)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_feature_plot_object <- function(input_rds) {
  readRDS(input_rds)
}

list_available_assays <- function(seurat_obj) {
  names(seurat_obj@assays)
}

list_available_reductions <- function(seurat_obj) {
  names(seurat_obj@reductions)
}

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

save_feature_plot <- function(plot_obj,
                              output_file,
                              width = 12,
                              height = 10) {
  ensure_dir(dirname(output_file))
  pdf(output_file, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

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
