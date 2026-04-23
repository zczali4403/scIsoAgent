suppressPackageStartupMessages({
  library(data.table)
})

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

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

infer_sample_name <- function(barcode) {
  parts <- strsplit(barcode, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    return(NA_character_)
  }
  parts[[1]]
}

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
