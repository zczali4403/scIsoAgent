# ------------------------------------------------------------------------------
# Script: draw_features.r
# Purpose: One-shot CLI wrapper for drawing selected gene/isoform features on embeddings.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
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
          "Rscript scripts/draw_features.r",
          "--input_rds=/path/to/seurat_object.rds",
          "--features=Sox9,Neurog3,Fev,Rbp4",
          "--output_file=/path/to/output.pdf",
          "[--reduction=UMAP]",
          "[--assay=RNA]"
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

# Brief: Create the parent directory for an output file when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the output file parent directory exists; no analysis object is returned.
ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

# Brief: Resolve a requested reduction name in a Seurat object.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   reduction_arg: User-requested dimensional reduction name.
# Returns: Resolved reduction name value or matching records.
resolve_reduction_name <- function(seurat_obj, reduction_arg) {
  reduction_names <- names(seurat_obj@reductions)
  if (length(reduction_names) == 0) {
    stop("No reductions were found in the Seurat object")
  }

  if (reduction_arg %in% reduction_names) {
    return(reduction_arg)
  }

  lower_map <- stats::setNames(reduction_names, tolower(reduction_names))
  lookup_key <- tolower(reduction_arg)
  if (lookup_key %in% names(lower_map)) {
    return(lower_map[[lookup_key]])
  }

  stop(
    paste(
      "Requested reduction was not found:",
      reduction_arg,
      "\nAvailable reductions:",
      paste(reduction_names, collapse = ", ")
    )
  )
}

# Brief: Stop if any requested feature is absent from the active assay.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
# Returns: No return value; raises an error if validation fails.
validate_features <- function(seurat_obj, features) {
  available_features <- rownames(seurat_obj)
  missing_features <- features[!features %in% available_features]

  if (length(missing_features) > 0) {
    stop(
      paste(
        "The following features were not found in the active assay:",
        paste(missing_features, collapse = ", ")
      )
    )
  }
}

# Brief: Sanitize text so it is safe for filenames or labels.
# Args:
#   x: Raw label text to sanitize.
# Returns: Constructed safe label object/table.
make_safe_label <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  if (x == "") {
    x <- "feature_plot"
  }
  x
}

# Brief: Write a compact run summary table for the current script.
# Args:
#   output_file: Destination path for the generated table or figure.
#   input_rds: Path to the Seurat RDS object used as input.
#   assay_name: Seurat assay name to set or validate.
#   reduction_name: Resolved dimensional reduction name.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
# Returns: Writes the module-specific summary TSV next to the main output file.
write_summary <- function(output_file, input_rds, assay_name, reduction_name, features) {
  summary_df <- data.frame(
    metric = c("input_rds", "assay", "reduction", "n_features", "features"),
    value = c(input_rds, assay_name, reduction_name, length(features), paste(features, collapse = ",")),
    stringsAsFactors = FALSE
  )
  fwrite(
    summary_df,
    file = sub("\\.pdf$", "_summary.tsv", output_file, ignore.case = TRUE),
    sep = "\t"
  )
}

# Brief: Draw and save a one-shot FeatureDimPlot for selected features.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   input_rds: Path to the Seurat RDS object used as input.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
#   output_file: Destination path for the generated table or figure.
#   reduction_arg: User-requested dimensional reduction name.
#   assay_name: Seurat assay name to set or validate.
# Returns: Writes the FeatureDimPlot PDF and companion summary TSV for the requested features.
draw_feature_plot <- function(seurat_obj,
                              input_rds,
                              features,
                              output_file,
                              reduction_arg = "UMAP",
                              assay_name = NULL) {
  working_obj <- seurat_obj

  if (!is.null(assay_name) && assay_name != "") {
    if (!assay_name %in% names(working_obj@assays)) {
      stop(
        paste(
          "Requested assay was not found:",
          assay_name,
          "\nAvailable assays:",
          paste(names(working_obj@assays), collapse = ", ")
        )
      )
    }
    DefaultAssay(working_obj) <- assay_name
  }

  reduction_name <- resolve_reduction_name(working_obj, reduction_arg)
  validate_features(working_obj, features)

  p <- FeatureDimPlot(
    srt = working_obj,
    features = features,
    reduction = reduction_name,
    theme_use = "theme_blank"
  ) +
    ggtitle(paste0("FeatureDimPlot on ", reduction_name))

  ensure_parent_dir(output_file)
  pdf(output_file, width = max(7, 3.5 * min(length(features), 4)), height = 6, bg = "white")
  print(p)
  dev.off()

  write_summary(
    output_file = output_file,
    input_rds = input_rds,
    assay_name = DefaultAssay(working_obj),
    reduction_name = reduction_name,
    features = features
  )
}

# Main workflow entrypoint.
# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_args(args, c("input_rds", "features", "output_file"))

  input_rds <- args$input_rds
  output_file <- args$output_file
  features <- split_csv_arg(args$features)
  reduction_arg <- if ("reduction" %in% names(args)) args$reduction else "UMAP"
  assay_name <- if ("assay" %in% names(args)) args$assay else NULL

  seurat_obj <- readRDS(input_rds)

  if (length(features) == 0) {
    stop("--features must contain at least one feature")
  }

  draw_feature_plot(
    seurat_obj = seurat_obj,
    input_rds = input_rds,
    features = features,
    output_file = output_file,
    reduction_arg = reduction_arg,
    assay_name = assay_name
  )
}

main()
