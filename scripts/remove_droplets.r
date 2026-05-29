# ------------------------------------------------------------------------------
# Script: remove_droplets.r
# Purpose: Run empty-droplet detection and generate filtered Seurat outputs.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SCP)
  library(DropletUtils)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(BiocParallel)
  library(data.table)
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
        paste(missing_keys, collapse = ", "),
        "\nRequired usage:",
        paste(
          "Rscript scripts/remove_droplets.r",
          "--gene_count_file=/path/to/gene_count.csv",
          "--output_path=/path/to/output_dir",
          "--output_seurat_file=sample_name",
          "[--fdr_threshold=0.001]",
          "[--lower=500]",
          "[--project_name=remove_droplets]"
        )
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

# Brief: Convert a dense count matrix into a sparse Seurat-compatible count matrix.
# Args:
#   count_matrix: Gene-by-cell count matrix used to create a sparse Seurat counts layer.
# Returns: Constructed dgcmatrix object/table.
make_dgcmatrix <- function(count_matrix) {
  seurat_object <- CreateSeuratObject(counts = count_matrix, project = "singlecell")
  GetAssayData(seurat_object, assay = "RNA", slot = "counts")
}

# Brief: Extract cell IDs from a clustered Seurat object for overlap checks with emptyDrops results.
# Args:
#   seurat_object: Seurat object being processed or plotted.
# Returns: Data frame containing cell_id values from Seurat cluster metadata.
overlap_true_cell <- function(seurat_object) {
  seurat_cluster.df <- as.data.frame(seurat_object$seurat_clusters)
  tibble::rownames_to_column(seurat_cluster.df, "cell_id")
}

# Brief: Convert emptyDrops FDR values into cell-versus-empty-droplet metadata labels.
# Args:
#   e.out: emptyDrops result object containing per-barcode FDR values.
#   fdr_threshold: FDR cutoff used to classify emptyDrops cell calls.
# Returns: Metadata data frame containing FDR and Cells/Empty_drops classification for each barcode.
true_cells_metadata <- function(e.out, fdr_threshold) {
  cells <- as.data.frame(e.out@rownames)
  fdr <- as.data.frame(e.out$FDR)
  tf_cells <- cbind(cells, fdr)
  tf_cells <- data.frame(tf_cells[, -1], row.names = tf_cells[, 1])
  data.table::setnames(tf_cells, c("FDR"))

  tf_cells %>%
    mutate(
      is.cell = case_when(
        FDR < fdr_threshold ~ "Cells",
        TRUE ~ "Empty_drops"
      )
    )
}

# Brief: Create a barcode rank plot showing knee and inflection points.
# Args:
#   br.out: barcodeRanks result object used to draw knee/inflection diagnostics.
# Returns: ggplot/plot object for rankplot.
build_rankplot <- function(br.out) {
  knee_value <- metadata(br.out)$knee
  inflection_value <- metadata(br.out)$inflection

  ggplot(br.out, aes(x = rank, y = total)) +
    geom_point(color = "#3B82F6", alpha = 0.45, size = 0.7) +
    geom_line(aes(y = fitted), color = "#D9485F", linewidth = 0.8) +
    geom_hline(yintercept = knee_value, color = "#0F766E", linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = inflection_value, color = "#A16207", linetype = "dashed", linewidth = 0.7) +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_log10(labels = scales::label_number()) +
    annotate(
      "label",
      x = max(br.out$rank, na.rm = TRUE),
      y = knee_value,
      label = "Knee",
      hjust = 1.05,
      vjust = -0.2,
      size = 3,
      fill = "white",
      color = "#0F766E",
      label.size = 0.15
    ) +
    annotate(
      "label",
      x = max(br.out$rank, na.rm = TRUE),
      y = inflection_value,
      label = "Inflection",
      hjust = 1.05,
      vjust = 1.2,
      size = 3,
      fill = "white",
      color = "#A16207",
      label.size = 0.15
    ) +
    labs(
      title = "Barcode Rank Plot",
      subtitle = "EmptyDrops input barcode distribution with fitted curve",
      x = "Barcode rank",
      y = "Total counts"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E5E7EB"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Brief: Create a table grob summarizing empty-droplet filtering results.
# Args:
#   fdr_threshold: FDR cutoff used to classify emptyDrops cell calls.
#   lower: Lower total-count threshold passed to emptyDrops.
#   count_true_cells: Cell barcode or cell-level metadata used for filtering or grouping.
#   count_empty_drops: Number of barcodes classified as empty droplets.
# Returns: Constructed summary grob object/table.
build_summary_grob <- function(fdr_threshold, lower, count_true_cells, count_empty_drops) {
  summary_table <- data.frame(
    Metric = c("FDR threshold", "Lower count cutoff", "Detected cells", "Removed empty droplets"),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops)
  )

  tableGrob(
    summary_table,
    rows = NULL,
    theme = ttheme_minimal(
      base_size = 10,
      core = list(
        fg_params = list(hjust = 0.5, x = 0.5),
        bg_params = list(fill = rep("#FFFFFF", nrow(summary_table)), col = "#D1D5DB")
      ),
      colhead = list(
        fg_params = list(fontface = "bold", col = "white", hjust = 0.5, x = 0.5),
        bg_params = list(fill = "#111827", col = "#111827")
      )
    )
  )
}

# Brief: Run emptyDrops, generate QC plots, and save pre/post-filter Seurat objects.
# Args:
#   output_path: Directory path used for this step’s output files.
#   gene_count_file: Gene-by-cell count matrix CSV used by downstream steps.
#   output_seurat_file: Filesystem input or output location used by this step.
#   fdr_threshold: FDR cutoff used to classify emptyDrops cell calls.
#   lower: Lower total-count threshold passed to emptyDrops.
#   project_name: Project label stored in newly created Seurat objects.
# Returns: Writes emptyDrops plots, summary CSV, and with/removed-empty Seurat RDS files.
perform_empty_drops_analysis <- function(output_path,
                                         gene_count_file,
                                         output_seurat_file,
                                         fdr_threshold = 0.001,
                                         lower = 100,
                                         project_name = "remove_droplets") {
  ensure_dir(output_path)

  dt <- fread(gene_count_file)
  gene_names <- dt[[1]]
  dt[[1]] <- NULL
  df <- as.data.frame(dt)
  rownames(df) <- gene_names

  seurat_obj <- CreateSeuratObject(counts = df, min.features = 20)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  outs.ddcmatrix <- make_dgcmatrix(df)
  br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
  e.out <- emptyDrops(
    outs.ddcmatrix,
    lower = lower,
    niters = 10000,
    test.ambient = TRUE,
    BPPARAM = SerialParam()
  )

  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>%
    filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")

  overlap_CR <- overlap_true_cell(seurat_obj)
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)

  cells_CR <- true_cells_metadata(e.out, fdr_threshold)
  cells_CR <- cells_CR[colnames(seurat_obj), , drop = FALSE]
  seurat_obj <- AddMetaData(seurat_obj, metadata = cells_CR)

  rankplot <- build_rankplot(br.out)

  cell_counts <- as.data.frame(table(seurat_obj@meta.data$is.cell))
  count_true_cells <- ifelse(
    length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
    cell_counts$Freq[cell_counts$Var1 == "Cells"],
    0
  )
  count_empty_drops <- ifelse(
    length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
    cell_counts$Freq[cell_counts$Var1 == "Empty_drops"],
    0
  )

  summary_grob <- build_summary_grob(
    fdr_threshold = fdr_threshold,
    lower = lower,
    count_true_cells = count_true_cells,
    count_empty_drops = count_empty_drops
  )

  cell_dim_plot <- CellDimPlot(
    srt = seurat_obj,
    group.by = "is.cell",
    reduction = "umap",
    title = "Cells vs Empty Droplets"
  ) +
    scale_color_manual(values = c("Cells" = "#2563EB", "Empty_drops" = "#94A3B8")) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(plot.title = element_text(face = "bold"))

  ncount_plot <- FeatureDimPlot(
    srt = seurat_obj,
    features = "nCount_RNA",
    reduction = "umap",
    theme_use = "theme_blank"
  ) +
    ggtitle("UMAP by nCount_RNA")

  nfeature_plot <- FeatureDimPlot(
    srt = seurat_obj,
    features = "nFeature_RNA",
    reduction = "umap",
    theme_use = "theme_blank"
  ) +
    ggtitle("UMAP by nFeature_RNA")

  combined_plot <- grid.arrange(
    rankplot,
    cell_dim_plot,
    ncount_plot,
    nfeature_plot,
    summary_grob,
    ncol = 2,
    layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)),
    top = textGrob(
      "Empty Droplet Removal Summary",
      gp = gpar(fontsize = 14, fontface = "bold")
    )
  )

  pdf(
    file = file.path(output_path, paste0(output_seurat_file, "_plots.pdf")),
    width = 11,
    height = 14,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(combined_plot)
  dev.off()

  summary_table <- data.frame(
    metric = c("fdr_threshold", "lower_count_cutoff", "detected_cells", "removed_empty_drops"),
    value = c(fdr_threshold, lower, count_true_cells, count_empty_drops)
  )
  write.csv(
    summary_table,
    file = file.path(output_path, paste0(output_seurat_file, "_summary.csv")),
    row.names = FALSE
  )

  seurat_obj_rm_empty <- subset(seurat_obj, subset = is.cell == "Cells")

  saveRDS(seurat_obj, file = file.path(output_path, paste0("with_empty_", output_seurat_file, ".rds")))
  saveRDS(seurat_obj_rm_empty, file = file.path(output_path, paste0("removed_empty_", output_seurat_file, ".rds")))
}

# Main workflow entrypoint.
# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required_keys <- c("gene_count_file", "output_path", "output_seurat_file")
  require_args(args, required_keys)

  fdr_threshold <- if ("fdr_threshold" %in% names(args)) as.numeric(args$fdr_threshold) else 0.001
  lower <- if ("lower" %in% names(args)) as.numeric(args$lower) else 500
  project_name <- if ("project_name" %in% names(args)) args$project_name else "remove_droplets"

  perform_empty_drops_analysis(
    output_path = args$output_path,
    gene_count_file = args$gene_count_file,
    output_seurat_file = args$output_seurat_file,
    fdr_threshold = fdr_threshold,
    lower = lower,
    project_name = project_name
  )
}

main()
