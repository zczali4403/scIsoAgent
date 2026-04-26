library(rtracklayer)
library(data.table)
library(tidyverse)
library(Matrix)

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

split_csv_arg <- function(value) {
  if (is.null(value) || value == "") {
    return(character(0))
  }

  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

require_common_args <- function(parsed) {
  common_keys <- c(
    "tool",
    "reference_gtf_file",
    "dict_output_file",
    "sample_names",
    "barcodes_paths",
    "features_paths",
    "matrix_paths",
    "matrix_output_file"
  )

  missing_keys <- common_keys[!common_keys %in% names(parsed)]
  if (length(missing_keys) > 0) {
    stop(
      paste(
        "Missing required arguments:",
        paste(missing_keys, collapse = ", ")
      )
    )
  }
}

normalize_sample_names <- function(sample_names) {
  normalized_sample_names <- gsub("_", "-", sample_names, fixed = TRUE)

  if (any(normalized_sample_names == "")) {
    stop("Failed to infer sample names from --barcodes_paths")
  }

  if (!identical(sample_names, normalized_sample_names)) {
    message(
      paste(
        "Normalized sample_names by replacing underscores with hyphens:",
        paste(
          paste0(sample_names, " -> ", normalized_sample_names),
          collapse = ", "
        )
      )
    )
  }

  if (anyDuplicated(normalized_sample_names)) {
    stop(
      paste(
        "Normalized sample names are not unique after replacing underscores with hyphens.",
        "Please provide unique --sample_names explicitly."
      )
    )
  }

  normalized_sample_names
}

validate_tool_args <- function(parsed) {
  tool <- parsed$tool

  if (!tool %in% c("isoquant", "flair")) {
    stop("--tool must be one of: isoquant, flair")
  }

  if (tool == "isoquant" && !"output_gtf_file" %in% names(parsed)) {
    stop("IsoQuant mode requires --output_gtf_file")
  }

  if (tool == "flair" && !"flair_counts_tsv" %in% names(parsed)) {
    stop("FLAIR mode requires --flair_counts_tsv")
  }
}

validate_sample_lists <- function(sample_names, barcodes_paths, features_paths, matrix_paths) {
  n <- length(sample_names)

  if (n == 0) {
    stop("At least one sample must be provided in --sample_names")
  }

  if (length(barcodes_paths) != n || length(features_paths) != n || length(matrix_paths) != n) {
    stop(
      paste(
        "sample_names, barcodes_paths, features_paths, and matrix_paths must have the same number of entries",
        "\nCounts:",
        paste(
          c(
            paste0("sample_names=", length(sample_names)),
            paste0("barcodes_paths=", length(barcodes_paths)),
            paste0("features_paths=", length(features_paths)),
            paste0("matrix_paths=", length(matrix_paths))
          ),
          collapse = ", "
        )
      )
    )
  }
}

make_isoquant_isoform_gene_symbol_dict <- function(output_gtf_file,
                                                   reference_gtf_file,
                                                   dict_output_file) {
  gtf1 <- import(output_gtf_file)
  gtf1_df <- as.data.frame(gtf1)
  gtf1_transcript <- gtf1_df[gtf1_df$type == "transcript", ]
  selected_columns1 <- unique(gtf1_transcript[, c("transcript_id", "gene_id")])

  gtf2 <- import(reference_gtf_file)
  gtf2_df <- as.data.frame(gtf2)
  selected_columns2 <- unique(gtf2_df[, c("gene_name", "gene_id")])

  combined_data <- merge(
    selected_columns1,
    selected_columns2,
    by = "gene_id",
    all.x = TRUE
  )

  combined_data$gene_symbol <- ifelse(
    is.na(combined_data$gene_name),
    combined_data$gene_id,
    combined_data$gene_name
  )

  final_combined_data <- combined_data[, c("transcript_id", "gene_id", "gene_symbol")]

  ensure_parent_dir(dict_output_file)
  write.csv(final_combined_data, file = dict_output_file, row.names = FALSE)

  final_combined_data
}

make_flair_isoform_gene_symbol_dict <- function(flair_counts_tsv,
                                                reference_gtf_file,
                                                dict_output_file) {
  flair_counts <- fread(flair_counts_tsv, header = TRUE)

  if (ncol(flair_counts) < 1) {
    stop("FLAIR counts file must contain at least one column")
  }

  transcript_ids <- flair_counts[[1]]
  transcript_gene_df <- data.frame(
    transcript_id = transcript_ids,
    gene_id = sub("^.*_", "", transcript_ids),
    stringsAsFactors = FALSE
  ) %>%
    distinct(transcript_id, .keep_all = TRUE)

  gtf2 <- import(reference_gtf_file)
  gtf2_df <- as.data.frame(gtf2)
  selected_columns2 <- unique(gtf2_df[, c("gene_name", "gene_id")])

  combined_data <- merge(
    transcript_gene_df,
    selected_columns2,
    by = "gene_id",
    all.x = TRUE
  )

  combined_data$gene_symbol <- ifelse(
    is.na(combined_data$gene_name),
    combined_data$gene_id,
    combined_data$gene_name
  )

  final_combined_data <- combined_data[, c("transcript_id", "gene_id", "gene_symbol")]

  ensure_parent_dir(dict_output_file)
  write.csv(final_combined_data, file = dict_output_file, row.names = FALSE)

  final_combined_data
}

make_isoform_gene_symbol_dict <- function(tool,
                                          reference_gtf_file,
                                          dict_output_file,
                                          output_gtf_file = NULL,
                                          flair_counts_tsv = NULL) {
  if (tool == "isoquant") {
    return(
      make_isoquant_isoform_gene_symbol_dict(
        output_gtf_file = output_gtf_file,
        reference_gtf_file = reference_gtf_file,
        dict_output_file = dict_output_file
      )
    )
  }

  if (tool == "flair") {
    return(
      make_flair_isoform_gene_symbol_dict(
        flair_counts_tsv = flair_counts_tsv,
        reference_gtf_file = reference_gtf_file,
        dict_output_file = dict_output_file
      )
    )
  }

  stop("Unsupported tool")
}

make_gene_count_matrix_from_mtx <- function(sample_name,
                                            barcodes_path,
                                            features_path,
                                            matrix_path) {
  barcodes <- fread(barcodes_path, header = FALSE)[[1]]
  prefixed_barcodes <- paste0(sample_name, "_", barcodes)

  features <- fread(features_path, header = FALSE)
  gene_ids <- features[[1]]

  mat <- readMM(matrix_path)

  if (nrow(mat) != length(gene_ids)) {
    stop(
      paste0(
        "matrix row count does not match features count for sample ",
        sample_name
      )
    )
  }

  if (ncol(mat) != length(barcodes)) {
    stop(
      paste0(
        "matrix column count does not match barcodes count for sample ",
        sample_name
      )
    )
  }

  mat_dense <- as.matrix(mat)
  gene_count_matrix <- data.frame(
    gene_id = gene_ids,
    mat_dense,
    check.names = FALSE
  )

  colnames(gene_count_matrix)[-1] <- prefixed_barcodes
  gene_count_matrix
}

merge_sample_gene_count_matrices <- function(sample_names,
                                             barcodes_paths,
                                             features_paths,
                                             matrix_paths) {
  sample_matrices <- vector("list", length(sample_names))

  for (i in seq_along(sample_names)) {
    sample_matrices[[i]] <- make_gene_count_matrix_from_mtx(
      sample_name = sample_names[[i]],
      barcodes_path = barcodes_paths[[i]],
      features_path = features_paths[[i]],
      matrix_path = matrix_paths[[i]]
    )
  }

  merged_matrix <- Reduce(
    function(x, y) merge(x, y, by = "gene_id", all = TRUE),
    sample_matrices
  )

  merged_matrix[is.na(merged_matrix)] <- 0
  merged_matrix
}

convert_ENSGID_to_geneSymbol <- function(gene_count_matrix,
                                         id_symbol_df,
                                         matrix_output_file,
                                         return_df = FALSE) {
  id_symbol_df <- as_tibble(id_symbol_df) %>%
    dplyr::select(gene_id, gene_symbol) %>%
    distinct(gene_id, .keep_all = TRUE)

  colnames(gene_count_matrix)[1] <- "gene_id"

  formatted_gene_count_matrix <- gene_count_matrix %>%
    merge(id_symbol_df, by = "gene_id") %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::select(-gene_id) %>%
    column_to_rownames(var = "gene_symbol")

  ensure_parent_dir(matrix_output_file)
  fwrite(formatted_gene_count_matrix, matrix_output_file, row.names = TRUE)

  if (return_df) {
    return(formatted_gene_count_matrix)
  }

  invisible(NULL)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_common_args(args)
  validate_tool_args(args)

  barcodes_paths <- split_csv_arg(args$barcodes_paths)
  features_paths <- split_csv_arg(args$features_paths)
  matrix_paths <- split_csv_arg(args$matrix_paths)
  sample_names <- split_csv_arg(args$sample_names)
  sample_names <- normalize_sample_names(sample_names)

  validate_sample_lists(sample_names, barcodes_paths, features_paths, matrix_paths)

  isoform_gene_dict <- make_isoform_gene_symbol_dict(
    tool = args$tool,
    reference_gtf_file = args$reference_gtf_file,
    dict_output_file = args$dict_output_file,
    output_gtf_file = args$output_gtf_file,
    flair_counts_tsv = args$flair_counts_tsv
  )

  merged_gene_count_matrix <- merge_sample_gene_count_matrices(
    sample_names = sample_names,
    barcodes_paths = barcodes_paths,
    features_paths = features_paths,
    matrix_paths = matrix_paths
  )

  convert_ENSGID_to_geneSymbol(
    gene_count_matrix = merged_gene_count_matrix,
    id_symbol_df = isoform_gene_dict,
    matrix_output_file = args$matrix_output_file
  )
}

main()
