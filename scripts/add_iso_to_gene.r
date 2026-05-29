# ------------------------------------------------------------------------------
# Script: add_iso_to_gene.r
# Purpose: Merge isoform assay and reduction into a merged gene-level Seurat object.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
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
          "Rscript scripts/add_iso_to_gene.r",
          "--gene_rds=/path/to/merged_gene.rds",
          "--isoform_rds=/path/to/merged_isoform.rds",
          "--output_rds=/path/to/output_with_isoform.rds",
          "[--isoform_assay_name=isoform]",
          "[--isoform_reduction_name=isoform_umap]"
        )
      )
    )
  }
}

# Brief: Create the parent directory for an output file when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the output file parent directory exists; no analysis object is returned.
ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

# Brief: Copy isoform UMAP embeddings into a reduction object aligned to shared cells.
# Args:
#   isoform_obj: Isoform-level Seurat object.
#   common_cells: Cell barcode or cell-level metadata used for filtering or grouping.
#   reduction_name: Resolved dimensional reduction name.
#   new_reduction_name: Reduction name used when copying isoform embeddings into another object.
# Returns: Constructed isoform reduction object/table.
build_isoform_reduction <- function(isoform_obj,
                                    common_cells,
                                    reduction_name = "umap",
                                    new_reduction_name = "isoform_umap") {
  if (!reduction_name %in% names(isoform_obj@reductions)) {
    return(NULL)
  }

  reduction_obj <- isoform_obj[[reduction_name]]
  embeddings <- Embeddings(reduction_obj)
  embeddings <- embeddings[common_cells, , drop = FALSE]

  global_flag <- tryCatch(reduction_obj@global, error = function(e) FALSE)

  CreateDimReducObject(
    embeddings = embeddings,
    key = Key(reduction_obj),
    assay = DefaultAssay(isoform_obj),
    global = global_flag
  )
}

# Brief: Write a compact run summary table for the current script.
# Args:
#   output_rds: Path where the output Seurat RDS object should be saved.
#   gene_rds: Merged gene-level Seurat RDS path.
#   isoform_rds: Merged isoform-level Seurat RDS path.
#   common_cells: Cell barcode or cell-level metadata used for filtering or grouping.
#   gene_obj: Gene-level Seurat object.
#   isoform_obj: Isoform-level Seurat object.
#   assay_name: Seurat assay name to set or validate.
#   reduction_added: Dimensional reduction name used for embedding operations.
# Returns: Writes the module-specific summary TSV next to the main output file.
write_summary <- function(output_rds,
                          gene_rds,
                          isoform_rds,
                          common_cells,
                          gene_obj,
                          isoform_obj,
                          assay_name,
                          reduction_added) {
  summary_df <- data.frame(
    metric = c(
      "gene_rds",
      "isoform_rds",
      "gene_cells_before_merge",
      "isoform_cells_before_merge",
      "shared_cells",
      "gene_features",
      "isoform_features",
      "isoform_assay_name",
      "isoform_reduction_added"
    ),
    value = c(
      gene_rds,
      isoform_rds,
      ncol(gene_obj),
      ncol(isoform_obj),
      length(common_cells),
      nrow(gene_obj),
      nrow(isoform_obj),
      assay_name,
      reduction_added
    ),
    stringsAsFactors = FALSE
  )

  fwrite(
    summary_df,
    file = sub("\\.rds$", "_summary.tsv", output_rds, ignore.case = TRUE),
    sep = "\t"
  )
}

# Main workflow entrypoint.
# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_args(args, c("gene_rds", "isoform_rds", "output_rds"))

  gene_rds <- args$gene_rds
  isoform_rds <- args$isoform_rds
  output_rds <- args$output_rds
  isoform_assay_name <- if ("isoform_assay_name" %in% names(args)) args$isoform_assay_name else "isoform"
  isoform_reduction_name <- if ("isoform_reduction_name" %in% names(args)) args$isoform_reduction_name else "isoform_umap"

  gene_obj <- readRDS(gene_rds)
  isoform_obj <- readRDS(isoform_rds)

  common_cells <- intersect(colnames(gene_obj), colnames(isoform_obj))
  if (length(common_cells) == 0) {
    stop("No overlapping cell barcodes were found between gene_rds and isoform_rds")
  }

  gene_obj <- subset(gene_obj, cells = common_cells)
  isoform_obj <- subset(isoform_obj, cells = common_cells)
  common_cells <- colnames(gene_obj)
  isoform_obj <- isoform_obj[, common_cells]

  isoform_counts <- GetAssayData(isoform_obj, assay = DefaultAssay(isoform_obj), slot = "counts")
  isoform_counts <- isoform_counts[, common_cells, drop = FALSE]

  gene_obj[[isoform_assay_name]] <- CreateAssayObject(counts = isoform_counts)

  isoform_reduction <- build_isoform_reduction(
    isoform_obj = isoform_obj,
    common_cells = common_cells,
    reduction_name = "umap",
    new_reduction_name = isoform_reduction_name
  )

  reduction_added <- "FALSE"
  if (!is.null(isoform_reduction)) {
    gene_obj[[isoform_reduction_name]] <- isoform_reduction
    reduction_added <- "TRUE"
  }

  ensure_parent_dir(output_rds)
  saveRDS(gene_obj, output_rds)

  write_summary(
    output_rds = output_rds,
    gene_rds = gene_rds,
    isoform_rds = isoform_rds,
    common_cells = common_cells,
    gene_obj = gene_obj,
    isoform_obj = isoform_obj,
    assay_name = isoform_assay_name,
    reduction_added = reduction_added
  )
}

main()
