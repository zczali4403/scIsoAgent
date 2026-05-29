# ------------------------------------------------------------------------------
# Script: split_gene_count_by_sample.r
# Purpose: Split merged gene-count matrices into per-sample files by barcode prefix.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
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
          "Rscript scripts/split_gene_count_by_sample.r",
          "--gene_count_file=/path/to/geneSymbol_gene_count.csv",
          "--output_dir=/path/to/output_dir",
          "[--strip_prefix=false]"
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

# Brief: Infer a sample name from a prefixed barcode string.
# Args:
#   barcode: Cell barcode or cell-level metadata used for filtering or grouping.
# Returns: Inferred sample name values.
infer_sample_name <- function(barcode) {
  parts <- strsplit(barcode, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    return(NA_character_)
  }
  parts[[1]]
}

# Brief: Split a merged gene-count matrix into one file per sample using barcode prefixes.
# Args:
#   gene_count_file: Gene-by-cell count matrix CSV used by downstream steps.
#   output_dir: Directory where output files should be written.
#   strip_prefix: Whether to remove sample prefixes from barcodes after splitting or extraction.
# Returns: Writes per-sample gene-count CSV files and split_manifest.tsv; returns invisibly.
split_gene_count_by_sample <- function(gene_count_file, output_dir, strip_prefix = FALSE) {
  ensure_dir(output_dir)

  dt <- fread(gene_count_file, check.names = FALSE)
  if (ncol(dt) < 2) {
    stop("Input gene count matrix must have at least one feature column and one barcode column")
  }

  feature_col <- colnames(dt)[1]
  barcode_cols <- colnames(dt)[-1]
  sample_names <- vapply(barcode_cols, infer_sample_name, character(1))

  if (any(is.na(sample_names))) {
    bad_cols <- barcode_cols[is.na(sample_names)]
    stop(
      paste(
        "Failed to infer sample names from these barcode columns:",
        paste(bad_cols, collapse = ", "),
        "\nExpected barcode format: sample_barcode"
      )
    )
  }

  split_map <- split(barcode_cols, sample_names)
  manifest <- data.frame(
    sample_name = names(split_map),
    n_barcodes = vapply(split_map, length, integer(1)),
    output_file = character(length(split_map)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(split_map)) {
    sample_name <- names(split_map)[[i]]
    sample_cols <- split_map[[i]]
    sample_dt <- dt[, c(feature_col, sample_cols), with = FALSE]

    if (strip_prefix) {
      renamed_cols <- sub(paste0("^", sample_name, "_"), "", sample_cols)
      setnames(sample_dt, old = sample_cols, new = renamed_cols)
    }

    output_file <- file.path(output_dir, paste0(sample_name, "_geneSymbol_gene_count.csv"))
    fwrite(sample_dt, output_file)
    manifest$output_file[[i]] <- output_file
  }

  fwrite(manifest, file.path(output_dir, "split_manifest.tsv"), sep = "\t")
}

# Main workflow entrypoint.
# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_args(args, c("gene_count_file", "output_dir"))

  strip_prefix <- if ("strip_prefix" %in% names(args)) {
    tolower(args$strip_prefix) %in% c("true", "1", "yes")
  } else {
    FALSE
  }

  split_gene_count_by_sample(
    gene_count_file = args$gene_count_file,
    output_dir = args$output_dir,
    strip_prefix = strip_prefix
  )
}

main()
