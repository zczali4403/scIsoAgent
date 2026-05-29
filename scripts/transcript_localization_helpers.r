# ------------------------------------------------------------------------------
# Script: transcript_localization_helpers.r
# Purpose: Helper utilities for DeepLocRNA transcript localization parsing and plotting.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(tidyr)
})

# Brief: Create an output directory when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the target directory exists; no analysis object is returned.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Brief: Write a ggplot-compatible plot object to a PDF file.
# Args:
#   plot_obj: Plot object to save.
#   file_path: Filesystem path for an input or output file.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Writes the provided plot object to a PDF file.
save_plot_pdf <- function(plot_obj, file_path, width = 9, height = 6) {
  ensure_dir(dirname(file_path))
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Parse a DeepLocRNA prediction table into transcript scores and predicted classes.
# Args:
#   prediction_file: DeepLocRNA prediction output file.
# Returns: Loaded deeplocrna predictions object/table from disk or package resources.
read_deeplocrna_predictions <- function(prediction_file) {
  lines <- readLines(prediction_file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) < 2) {
    stop("DeepLocRNA prediction file appears malformed: expected a header plus at least one data row")
  }

  # DeepLocRNA output does not include an explicit first-column header for transcript ids.
  # Prepending one preserves the Nucleus column instead of accidentally renaming it away.
  lines[[1]] <- paste("transcript_id", lines[[1]])

  pred_df <- utils::read.table(
    text = lines,
    header = TRUE,
    sep = "",
    quote = '"',
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )

  if (ncol(pred_df) < 3) {
    stop("DeepLocRNA prediction file appears malformed: expected transcript ID, score columns, and Prediction")
  }

  pred_df$transcript_id <- gsub("^>", "", as.character(pred_df$transcript_id))

  required_cols <- c(
    "transcript_id",
    "Nucleus",
    "Exosome",
    "Cytosol",
    "Ribosome",
    "Membrane",
    "Endoplasmic reticulum",
    "Prediction"
  )
  missing_cols <- required_cols[!required_cols %in% colnames(pred_df)]
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "DeepLocRNA prediction file is missing required columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  score_cols <- setdiff(colnames(pred_df), c("transcript_id", "Prediction"))
  if (length(score_cols) == 0) {
    stop("No localization score columns were detected in the DeepLocRNA prediction file")
  }

  pred_df[score_cols] <- lapply(pred_df[score_cols], as.numeric)
  attr(pred_df, "score_cols") <- score_cols
  pred_df
}

# Brief: Return DeepLocRNA compartment score columns.
# Args:
#   pred_df: DeepLocRNA prediction data frame.
# Returns: Retrieved localization score columns values.
get_localization_score_columns <- function(pred_df) {
  setdiff(colnames(pred_df), c("transcript_id", "Prediction"))
}

# Brief: Count transcripts assigned to each DeepLocRNA prediction class.
# Args:
#   pred_df: DeepLocRNA prediction data frame.
# Returns: Summary table for prediction counts.
summarize_prediction_counts <- function(pred_df) {
  pred_df %>%
    dplyr::count(Prediction, name = "transcript_count", sort = TRUE)
}

# Brief: Plot transcript counts per DeepLocRNA localization class.
# Args:
#   summary_df: Summary table used for plotting or writing results.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for prediction class counts; writes a PDF when an output path is provided.
plot_prediction_class_counts <- function(summary_df, output_file = NULL) {
  plot_df <- summary_df %>%
    mutate(
      Prediction = reorder(Prediction, transcript_count),
      label_y = transcript_count + max(transcript_count, na.rm = TRUE) * 0.03
    )

  pred_palette <- c(
    "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
    "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#2F4B7C",
    "#D45087", "#FFA600"
  )
  palette_vals <- stats::setNames(
    pred_palette[seq_len(nrow(plot_df))],
    as.character(plot_df$Prediction)
  )

  p <- ggplot(plot_df, aes(x = Prediction, y = transcript_count, fill = Prediction)) +
    geom_col(
      color = "white",
      linewidth = 0.35,
      width = 0.78,
      show.legend = FALSE
    ) +
    geom_text(
      aes(y = label_y, label = transcript_count),
      hjust = 0,
      size = 3.4,
      color = "#243447"
    ) +
    scale_fill_manual(values = palette_vals) +
    coord_flip(clip = "off") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = "Transcript Counts by DeepLocRNA Prediction",
      x = "Prediction class",
      y = "Transcript count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2933"),
      axis.title = element_text(face = "bold", color = "#243447"),
      axis.text = element_text(color = "#243447"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#D9E2EC", linewidth = 0.35),
      plot.margin = margin(12, 24, 12, 12)
    )

  if (!is.null(output_file) && output_file != "") {
    plot_height <- max(6, 0.45 * nrow(plot_df) + 2.2)
    save_plot_pdf(p, output_file, width = 10.5, height = plot_height)
  }

  p
}

# Brief: Read the transcript-to-gene mapping table produced during resource setup.
# Args:
#   mapping_file: Transcript-to-gene mapping CSV produced by resource setup.
# Returns: Loaded transcript gene mapping object/table from disk or package resources.
read_transcript_gene_mapping <- function(mapping_file) {
  dict_df <- data.table::fread(mapping_file, data.table = FALSE)

  required_cols <- c("transcript_id", "gene_id", "gene_symbol")
  missing_cols <- required_cols[!required_cols %in% colnames(dict_df)]
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "The transcript-to-gene mapping file is missing required columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  dict_df %>%
    mutate(
      transcript_id = as.character(transcript_id),
      gene_id = as.character(gene_id),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    distinct(transcript_id, .keep_all = TRUE)
}

# Brief: Find all transcript records matching a gene symbol or gene ID.
# Args:
#   gene_query: Gene symbol or gene ID requested by the user.
#   mapping_df: Loaded transcript-to-gene mapping data frame.
# Returns: Resolved gene transcripts value or matching records.
resolve_gene_transcripts <- function(gene_query, mapping_df) {
  if (is.null(gene_query) || gene_query == "") {
    stop("gene_query must be a non-empty gene symbol or gene_id")
  }

  matches <- mapping_df %>%
    filter(
      gene_symbol == gene_query |
      gene_id == gene_query
    )

  if (nrow(matches) == 0) {
    stop(paste("No transcripts found for gene query:", gene_query))
  }

  matches %>%
    arrange(gene_symbol, gene_id, transcript_id)
}

# Brief: Write selected transcript IDs to a plain-text file.
# Args:
#   transcript_df: Transcript mapping subset for one selected gene.
#   output_file: Destination path for the generated table or figure.
# Returns: Writes transcript ID text file and invisibly returns the output path.
write_transcript_ids_file <- function(transcript_df, output_file) {
  if (!"transcript_id" %in% colnames(transcript_df)) {
    stop("transcript_df must contain a transcript_id column")
  }

  ensure_dir(dirname(output_file))
  writeLines(unique(as.character(transcript_df$transcript_id)), con = output_file)
  invisible(output_file)
}

# Brief: Build an AWK command that extracts selected transcript records from FASTA.
# Args:
#   transcript_ids_file: Plain-text file listing transcript IDs to extract from FASTA.
#   input_fasta: Input transcript FASTA file.
#   output_fasta: Output FASTA file to create.
# Returns: Shell command string ready to inspect and execute.
prepare_subset_fasta_command <- function(transcript_ids_file,
                                         input_fasta,
                                         output_fasta) {
  paste0(
    "awk '\n",
    "BEGIN { while ((getline < ", shQuote(transcript_ids_file), ") > 0) ids[$1] = 1 }\n",
    "/^>/ {\n",
    "  header = substr($0, 2)\n",
    "  split(header, parts, /[[:space:]]+/)\n",
    "  keep = (parts[1] in ids)\n",
    "}\n",
    "keep { print }\n",
    "' ", shQuote(input_fasta), " > ", shQuote(output_fasta)
  )
}

# Brief: Extract selected transcript FASTA records using an ID list.
# Args:
#   transcript_ids_file: Plain-text file listing transcript IDs to extract from FASTA.
#   input_fasta: Input transcript FASTA file.
#   output_fasta: Output FASTA file to create.
# Returns: List with output FASTA path and number of transcripts kept.
extract_subset_fasta_by_transcript_ids <- function(transcript_ids_file,
                                                   input_fasta,
                                                   output_fasta) {
  if (!file.exists(transcript_ids_file)) {
    stop(paste("Transcript ID file not found:", transcript_ids_file))
  }
  if (!file.exists(input_fasta)) {
    stop(paste("Input FASTA not found:", input_fasta))
  }

  ids <- readLines(transcript_ids_file, warn = FALSE)
  ids <- unique(trimws(ids[nzchar(trimws(ids))]))
  if (length(ids) == 0) {
    stop("No transcript IDs found in transcript_ids_file")
  }

  ensure_dir(dirname(output_fasta))

  keep <- FALSE
  n_kept <- 0L
  con_in <- file(input_fasta, open = "r")
  con_out <- file(output_fasta, open = "w")
  on.exit({
    try(close(con_in), silent = TRUE)
    try(close(con_out), silent = TRUE)
  }, add = TRUE)

  repeat {
    line <- readLines(con_in, n = 1, warn = FALSE)
    if (length(line) == 0) break

    if (startsWith(line, ">")) {
      header <- sub("^>", "", line)
      tx_id <- strsplit(header, "[[:space:]]+")[[1]][1]
      keep <- tx_id %in% ids
      if (keep) n_kept <- n_kept + 1L
    }

    if (keep) writeLines(line, con_out)
  }

  invisible(list(output_fasta = output_fasta, kept_transcripts = n_kept))
}

# Brief: Join DeepLocRNA predictions with mapping records for one selected gene.
# Args:
#   gene_query: Gene symbol or gene ID requested by the user.
#   pred_df: DeepLocRNA prediction data frame.
#   mapping_df: Loaded transcript-to-gene mapping data frame.
#   keep_prediction_col: Whether to keep the predicted DeepLocRNA class in gene-level output tables.
# Returns: List containing gene transcript table, score matrix, and score-column names.
prepare_gene_localization_scores <- function(gene_query,
                                             pred_df,
                                             mapping_df,
                                             keep_prediction_col = TRUE) {
  gene_tx <- resolve_gene_transcripts(gene_query, mapping_df)
  score_cols <- get_localization_score_columns(pred_df)

  merged_df <- gene_tx %>%
    inner_join(pred_df, by = "transcript_id")

  if (nrow(merged_df) == 0) {
    stop(paste("No DeepLocRNA predictions matched the transcripts for gene:", gene_query))
  }

  keep_cols <- c("transcript_id", "gene_id", "gene_symbol", score_cols)
  if (keep_prediction_col) {
    keep_cols <- c(keep_cols, "Prediction")
  }

  merged_df <- merged_df[, keep_cols, drop = FALSE]

  score_mat <- as.matrix(merged_df[, score_cols, drop = FALSE])
  rownames(score_mat) <- merged_df$transcript_id

  list(
    transcript_table = merged_df,
    score_matrix = score_mat,
    score_columns = score_cols
  )
}

# Brief: Write a gene-specific transcript localization score table.
# Args:
#   localization_obj: List containing gene-specific localization table and score matrix.
#   output_file: Destination path for the generated table or figure.
# Returns: Writes gene-level localization score CSV.
write_gene_localization_table <- function(localization_obj, output_file) {
  ensure_dir(dirname(output_file))
  write.csv(localization_obj$transcript_table, output_file, row.names = FALSE)
}

# Brief: Plot gene-specific DeepLocRNA scores relative to compartment thresholds.
# Args:
#   localization_obj: List containing gene-specific localization table and score matrix.
#   gene_query: Gene symbol or gene ID requested by the user.
#   output_file: Destination path for the generated table or figure.
#   cluster_rows: Whether heatmap rows should be clustered.
#   cluster_cols: Whether heatmap columns should be clustered.
# Returns: Plot object for gene localization heatmap; writes a PDF when an output path is provided.
plot_gene_localization_heatmap <- function(localization_obj,
                                           gene_query,
                                           output_file = NULL,
                                           cluster_rows = FALSE,
                                           cluster_cols = FALSE) {
  score_mat <- localization_obj$score_matrix

  annotation_row <- NULL
  annotation_colors <- NULL
  if ("Prediction" %in% colnames(localization_obj$transcript_table)) {
    annotation_row <- data.frame(
      Prediction = localization_obj$transcript_table$Prediction,
      row.names = localization_obj$transcript_table$transcript_id,
      stringsAsFactors = FALSE
    )

    pred_levels <- unique(annotation_row$Prediction)
    pred_palette <- grDevices::colorRampPalette(
      c("#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756", "#72B7B2")
    )(length(pred_levels))
    annotation_colors <- list(
      Prediction = stats::setNames(pred_palette, pred_levels)
    )
  }

  threshold_map <- c(
    "Nucleus" = 0.7551,
    "Exosome" = 0.9796,
    "Cytosol" = 0.2245,
    "Ribosome" = 0.2857,
    "Membrane" = 0.3061,
    "Endoplasmic reticulum" = 0.1837
  )

  keep_cols <- intersect(names(threshold_map), colnames(score_mat))
  if (length(keep_cols) == 0) {
    stop("None of the threshold-defined DeepLocRNA compartments were found in the score matrix")
  }

  score_mat <- score_mat[, keep_cols, drop = FALSE]
  col_thresholds <- threshold_map[keep_cols]
  col_labels <- paste0(
    keep_cols,
    "
(thr=",
    sprintf("%.3f", col_thresholds),
    ")"
  )
  colnames(score_mat) <- col_labels

  relative_margin_mat <- matrix(
    NA_real_,
    nrow = nrow(score_mat),
    ncol = ncol(score_mat),
    dimnames = dimnames(score_mat)
  )
  display_numbers <- matrix(
    "",
    nrow = nrow(score_mat),
    ncol = ncol(score_mat),
    dimnames = dimnames(score_mat)
  )

  for (i in seq_along(col_thresholds)) {
    thr <- col_thresholds[[i]]
    if (!is.na(thr) && thr > 0) {
      margin_vals <- score_mat[, i] / thr - 1
      relative_margin_mat[, i] <- margin_vals
      display_numbers[, i] <- ifelse(margin_vals > 0, "PASS", "")
    }
  }

  finite_vals <- relative_margin_mat[is.finite(relative_margin_mat)]
  min_margin <- if (length(finite_vals) > 0) {
    min(finite_vals, na.rm = TRUE)
  } else {
    -0.25
  }
  max_margin <- if (length(finite_vals) > 0) {
    max(finite_vals, na.rm = TRUE)
  } else {
    0.25
  }
  min_margin <- min(min_margin, -0.02)
  max_margin <- max(max_margin, 0.02)

  neg_palette <- grDevices::colorRampPalette(
    c("#24557A", "#7DB7D8", "#EAF4FF")
  )(80)
  pos_palette <- grDevices::colorRampPalette(
    c("#F6D38B", "#E58F65", "#B23A48")
  )(80)
  margin_palette <- c(neg_palette, pos_palette)
  neg_breaks <- seq(
    min_margin,
    0,
    length.out = length(neg_palette) + 1
  )
  pos_breaks <- seq(
    0,
    max_margin,
    length.out = length(pos_palette) + 1
  )
  breaks_vec <- c(neg_breaks, pos_breaks[-1])

  plot_width <- max(11.5, 1.85 * ncol(score_mat) + 3.8)
  plot_height <- max(6.5, 0.72 * nrow(score_mat) + 3.2)
  heat_cellwidth <- max(34, min(58, 250 / max(1, ncol(score_mat))))
  heat_cellheight <- max(20, min(34, 300 / max(1, nrow(score_mat))))
  legend_breaks <- c(min_margin, 0, max_margin)
  legend_labels <- c("Below", "At threshold", "Above")

  # Brief: Draw the already prepared DeepLocRNA heatmap inside the enclosing function.
  # Args: None.
  # Returns: Renders the nested DeepLocRNA pheatmap for the prepared localization matrix.
  draw_heatmap <- function() {
    pheatmap::pheatmap(
      relative_margin_mat,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      color = margin_palette,
      breaks = breaks_vec,
      border_color = "#FFFFFF",
      main = paste0(
        "Transcript Localization Relative To Threshold: ",
        gene_query,
        "\nNegative = below threshold, positive = above threshold"
      ),
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      display_numbers = display_numbers,
      number_color = "#202124",
      fontsize_number = 8.5,
      fontsize_row = 9.5,
      fontsize_col = 10.5,
      angle_col = "45",
      cellwidth = heat_cellwidth,
      cellheight = heat_cellheight,
      treeheight_row = if (cluster_rows) 42 else 0,
      treeheight_col = if (cluster_cols) 36 else 0,
      legend = TRUE,
      legend_breaks = legend_breaks,
      legend_labels = legend_labels,
      fontsize = 10.5,
      na_col = "#F3F6FA"
    )
  }

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    pdf(output_file, width = plot_width, height = plot_height, bg = "white")
    draw_heatmap()
    dev.off()
    return(invisible(NULL))
  }

  draw_heatmap()
}
