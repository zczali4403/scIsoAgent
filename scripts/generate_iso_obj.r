# ------------------------------------------------------------------------------
# Script: generate_iso_obj.r
# Purpose: Build per-sample isoform-level Seurat objects and subset by QC-retained cells.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
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
          "Rscript scripts/generate_iso_obj.r",
          "--barcodes_paths=/path/s1_barcodes.tsv,/path/s2_barcodes.tsv",
          "--features_paths=/path/s1_features.tsv,/path/s2_features.tsv",
          "--matrix_paths=/path/s1_matrix.mtx,/path/s2_matrix.mtx",
          "--gene_rds_paths=/path/s1_gene_singlets.rds,/path/s2_gene_singlets.rds",
          "--output_dir=/path/to/output_dir",
          "[--sample_names=sample1,sample2]",
          "[--min_features=0]"
        )
      )
    )
  }
}

# Brief: Convert a comma-separated command-line value into a trimmed character vector.
# Args:
#   value: Input string or scalar value to parse or summarize.
# Returns: Character vector/list of trimmed values from the comma-separated input.
split_csv_arg <- function(value) {
  if (is.null(value) || value == "") {
    return(character(0))
  }

  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

# Brief: Create an output directory when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the target directory exists; no analysis object is returned.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Brief: Infer sample names from input filenames when explicit sample labels are not provided.
# Args:
#   barcodes_paths: Per-sample barcode TSV paths, aligned with sample names.
# Returns: Inferred sample names values.
infer_sample_names <- function(barcodes_paths) {
  sample_names <- basename(barcodes_paths)
  sample_names <- sub("\\.barcodes\\.tsv(\\.gz)?$", "", sample_names)
  sample_names <- sub("\\.tsv(\\.gz)?$", "", sample_names, ignore.case = TRUE)
  sample_names <- sub("\\.txt(\\.gz)?$", "", sample_names, ignore.case = TRUE)
  sample_names <- gsub("[^A-Za-z0-9._-]", "_", sample_names)

  if (any(sample_names == "")) {
    stop("Failed to infer sample names from --barcodes_paths")
  }

  if (anyDuplicated(sample_names)) {
    stop("Inferred sample names are not unique; please provide --sample_names explicitly")
  }

  sample_names
}

# Brief: Verify that sample names and per-sample input file lists have matching lengths.
# Args:
#   sample_names: Sample labels aligned with the per-sample inputs.
#   barcodes_paths: Per-sample barcode TSV paths, aligned with sample names.
#   features_paths: Per-sample feature TSV paths, aligned with sample names.
#   matrix_paths: Per-sample MatrixMarket MTX paths, aligned with sample names.
#   gene_rds_paths: Per-sample gene-level Seurat RDS paths used for isoform filtering.
# Returns: No return value; raises an error if validation fails.
validate_sample_lists <- function(sample_names,
                                  barcodes_paths,
                                  features_paths,
                                  matrix_paths,
                                  gene_rds_paths) {
  n <- length(sample_names)

  if (n == 0) {
    stop("At least one sample must be provided")
  }

  lengths_to_check <- c(
    barcodes_paths = length(barcodes_paths),
    features_paths = length(features_paths),
    matrix_paths = length(matrix_paths),
    gene_rds_paths = length(gene_rds_paths)
  )

  if (any(lengths_to_check != n)) {
    stop(
      paste(
        "sample_names, barcodes_paths, features_paths, matrix_paths, and gene_rds_paths must have the same number of entries",
        "\nCounts:",
        paste(
          c(paste0("sample_names=", n), paste0(names(lengths_to_check), "=", lengths_to_check)),
          collapse = ", "
        )
      )
    )
  }
}

# Brief: Read feature identifiers from a matrix feature file.
# Args:
#   features_path: Feature TSV path for one matrix-market triplet.
# Returns: Loaded feature ids object/table from disk or package resources.
read_feature_ids <- function(features_path) {
  features_dt <- fread(features_path, header = FALSE)
  if (ncol(features_dt) < 1) {
    stop(paste("Feature file is empty:", features_path))
  }
  as.character(features_dt[[1]])
}

# Brief: Read cell barcode identifiers from a matrix barcode file.
# Args:
#   barcodes_path: Barcode TSV path for one matrix-market triplet.
# Returns: Loaded barcodes object/table from disk or package resources.
read_barcodes <- function(barcodes_path) {
  barcodes_dt <- fread(barcodes_path, header = FALSE)
  if (ncol(barcodes_dt) < 1) {
    stop(paste("Barcode file is empty:", barcodes_path))
  }
  as.character(barcodes_dt[[1]])
}

# Brief: Create one isoform-level Seurat object from barcode, feature, and matrix files.
# Args:
#   sample_name: Sample label used for barcode prefixing and output naming.
#   barcodes_path: Barcode TSV path for one matrix-market triplet.
#   features_path: Feature TSV path for one matrix-market triplet.
#   matrix_path: MatrixMarket MTX path for one count matrix.
#   min_features: Minimum detected features required to keep a cell or create an object.
# Returns: Isoform-level Seurat object for one sample with sample-prefixed barcodes.
create_isoform_seurat_object <- function(sample_name,
                                         barcodes_path,
                                         features_path,
                                         matrix_path,
                                         min_features = 0) {
  barcodes <- read_barcodes(barcodes_path)
  features <- read_feature_ids(features_path)
  mat <- readMM(matrix_path)

  if (nrow(mat) != length(features)) {
    stop(paste0("matrix row count does not match feature count for sample ", sample_name))
  }

  if (ncol(mat) != length(barcodes)) {
    stop(paste0("matrix column count does not match barcode count for sample ", sample_name))
  }

  prefixed_barcodes <- paste0(sample_name, "_", barcodes)
  rownames(mat) <- features
  colnames(mat) <- prefixed_barcodes

  CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.features = min_features
  )
}

# Brief: Create violin plots for QC metrics or isoform-level count summaries.
# Args:
#   seurat_object: Seurat object being processed or plotted.
#   sample_name: Sample label used for barcode prefixing and output naming.
# Returns: ggplot/plot object for violin plot.
build_violin_plot <- function(seurat_object, sample_name) {
  md <- seurat_object@meta.data[, c("nFeature_RNA", "nCount_RNA"), drop = FALSE]
  long_df <- stack(md)
  colnames(long_df) <- c("value", "metric")

  ggplot(long_df, aes(x = metric, y = value, fill = metric)) +
    geom_violin(scale = "width", trim = TRUE, color = "#374151", linewidth = 0.2) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    labs(
      title = paste0(sample_name, " isoform QC metrics before subsetting"),
      x = NULL,
      y = "Value"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# Brief: Save an isoform QC violin plot as a PDF file.
# Args:
#   plot_obj: Plot object to save.
#   file_path: Filesystem path for an input or output file.
# Returns: Writes the isoform QC violin plot PDF.
save_violin_plot <- function(plot_obj, file_path) {
  pdf(file_path, width = 10, height = 4.8, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Restrict an isoform Seurat object to cells retained in the matched gene-level QC object.
# Args:
#   isoform_obj: Isoform-level Seurat object.
#   gene_rds_path: Gene-level Seurat RDS used to define retained cells for one sample.
#   sample_name: Sample label used for barcode prefixing and output naming.
# Returns: Isoform Seurat object subset to cells shared with the matched gene-level object.
subset_isoform_object_by_gene_object <- function(isoform_obj, gene_rds_path, sample_name) {
  gene_obj <- readRDS(gene_rds_path)
  keep_cells <- intersect(colnames(isoform_obj), colnames(gene_obj))

  if (length(keep_cells) == 0) {
    stop(
      paste(
        "No overlapping cell barcodes were found between the isoform object and gene Seurat object for sample",
        sample_name
      )
    )
  }

  subset(isoform_obj, cells = keep_cells)
}

# Brief: Write per-sample input, cell-count, and feature-count summary tables.
# Args:
#   summary_rows: Per-sample summary records collected during isoform object generation.
#   output_dir: Directory where output files should be written.
# Returns: Writes a per-sample summary TSV for the current module.
write_sample_summary <- function(summary_rows, output_dir) {
  summary_df <- do.call(rbind, lapply(summary_rows, as.data.frame, stringsAsFactors = FALSE))
  fwrite(summary_df, file.path(output_dir, "generate_iso_obj_summary.tsv"), sep = "\t")
}

# Main workflow entrypoint.
# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_args(
    args,
    c("barcodes_paths", "features_paths", "matrix_paths", "gene_rds_paths", "output_dir")
  )

  barcodes_paths <- split_csv_arg(args$barcodes_paths)
  features_paths <- split_csv_arg(args$features_paths)
  matrix_paths <- split_csv_arg(args$matrix_paths)
  gene_rds_paths <- split_csv_arg(args$gene_rds_paths)
  output_dir <- args$output_dir
  min_features <- if ("min_features" %in% names(args)) as.integer(args$min_features) else 0L

  sample_names <- if ("sample_names" %in% names(args)) {
    split_csv_arg(args$sample_names)
  } else {
    infer_sample_names(barcodes_paths)
  }

  validate_sample_lists(
    sample_names = sample_names,
    barcodes_paths = barcodes_paths,
    features_paths = features_paths,
    matrix_paths = matrix_paths,
    gene_rds_paths = gene_rds_paths
  )

  ensure_dir(output_dir)

  summary_rows <- vector("list", length(sample_names))

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[[i]]
    sample_output_dir <- file.path(output_dir, sample_name)
    ensure_dir(sample_output_dir)

    isoform_obj <- create_isoform_seurat_object(
      sample_name = sample_name,
      barcodes_path = barcodes_paths[[i]],
      features_path = features_paths[[i]],
      matrix_path = matrix_paths[[i]],
      min_features = min_features
    )

    violin_plot <- build_violin_plot(isoform_obj, sample_name)
    save_violin_plot(
      violin_plot,
      file.path(sample_output_dir, paste0(sample_name, "_isoform_violinplot.pdf"))
    )

    subset_iso_obj <- subset_isoform_object_by_gene_object(
      isoform_obj = isoform_obj,
      gene_rds_path = gene_rds_paths[[i]],
      sample_name = sample_name
    )

    saveRDS(
      subset_iso_obj,
      file.path(sample_output_dir, paste0(sample_name, "_isoform_subset.rds"))
    )

    summary_rows[[i]] <- list(
      sample_name = sample_name,
      isoform_cells_before_subset = ncol(isoform_obj),
      isoform_cells_after_subset = ncol(subset_iso_obj),
      isoform_features = nrow(subset_iso_obj),
      gene_rds_path = gene_rds_paths[[i]]
    )
  }

  write_sample_summary(summary_rows, output_dir)
}

main()
