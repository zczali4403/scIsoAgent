suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(tidyr)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_plot_pdf <- function(plot_obj, file_path, width = 9, height = 6) {
  ensure_dir(dirname(file_path))
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

read_deeplocrna_predictions <- function(prediction_file) {
  pred_df <- data.table::fread(prediction_file, data.table = FALSE)

  if (ncol(pred_df) < 3) {
    stop("DeepLocRNA prediction file appears malformed: expected transcript ID, score columns, and Prediction")
  }

  first_col <- colnames(pred_df)[1]
  colnames(pred_df)[1] <- "transcript_id"
  pred_df$transcript_id <- gsub("^>", "", as.character(pred_df$transcript_id))

  required_cols <- c("Prediction")
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
  attr(pred_df, "original_first_col") <- first_col
  pred_df
}

get_localization_score_columns <- function(pred_df) {
  setdiff(colnames(pred_df), c("transcript_id", "Prediction"))
}

summarize_prediction_counts <- function(pred_df) {
  pred_df %>%
    dplyr::count(Prediction, name = "transcript_count", sort = TRUE)
}

plot_prediction_class_counts <- function(summary_df, output_file = NULL) {
  plot_df <- summary_df %>%
    mutate(Prediction = reorder(Prediction, transcript_count))

  p <- ggplot(plot_df, aes(x = Prediction, y = transcript_count)) +
    geom_col(fill = "#4C78A8", color = "black", width = 0.8) +
    coord_flip() +
    labs(
      title = "Transcript Counts by DeepLocRNA Prediction",
      x = "Prediction class",
      y = "Transcript count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 9, height = 6)
  }

  p
}

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

write_transcript_ids_file <- function(transcript_df, output_file) {
  if (!"transcript_id" %in% colnames(transcript_df)) {
    stop("transcript_df must contain a transcript_id column")
  }

  ensure_dir(dirname(output_file))
  writeLines(unique(as.character(transcript_df$transcript_id)), con = output_file)
  invisible(output_file)
}

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

write_gene_localization_table <- function(localization_obj, output_file) {
  ensure_dir(dirname(output_file))
  write.csv(localization_obj$transcript_table, output_file, row.names = FALSE)
}

plot_gene_localization_heatmap <- function(localization_obj,
                                           gene_query,
                                           output_file = NULL,
                                           cluster_rows = FALSE,
                                           cluster_cols = FALSE) {
  score_mat <- localization_obj$score_matrix

  annotation_row <- NULL
  if ("Prediction" %in% colnames(localization_obj$transcript_table)) {
    annotation_row <- data.frame(
      Prediction = localization_obj$transcript_table$Prediction,
      row.names = localization_obj$transcript_table$transcript_id,
      stringsAsFactors = FALSE
    )
  }

  palette_fun <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))
  heat_colors <- palette_fun(100)

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    pdf(output_file, width = 10, height = max(4, 0.35 * nrow(score_mat) + 2), bg = "white")
    pheatmap::pheatmap(
      score_mat,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      color = heat_colors,
      border_color = "grey85",
      main = paste0("Transcript Localization Scores: ", gene_query),
      annotation_row = annotation_row,
      fontsize_row = 9,
      fontsize_col = 10,
      angle_col = 45
    )
    dev.off()
    return(invisible(NULL))
  }

  pheatmap::pheatmap(
    score_mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color = heat_colors,
    border_color = "grey85",
    main = paste0("Transcript Localization Scores: ", gene_query),
    annotation_row = annotation_row,
    fontsize_row = 9,
    fontsize_col = 10,
    angle_col = 45
  )
}
