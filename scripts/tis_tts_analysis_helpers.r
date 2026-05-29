# ------------------------------------------------------------------------------
# Script: tis_tts_analysis_helpers.r
# Purpose: Helper utilities for TranslationAI TIS/TTS parsing and visualization.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
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
save_plot_pdf <- function(plot_obj, file_path, width = 10, height = 6) {
  ensure_dir(dirname(file_path))
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Parse TranslationAI ORF output into transcript TIS/TTS coordinates and scores.
# Args:
#   orf_file: TranslationAI ORF prediction output file.
# Returns: Loaded translationai orf output object/table from disk or package resources.
read_translationai_orf_output <- function(orf_file) {
  lines <- readLines(orf_file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0) {
    stop("TranslationAI ORF output is empty")
  }

  parsed <- lapply(lines, function(line) {
    if (!startsWith(line, ">")) {
      stop("Unexpected TranslationAI ORF format: each non-empty line should start with '>'")
    }

    body <- substring(line, 2)
    parts <- strsplit(body, "\t")[[1]]
    if (length(parts) < 2) {
      stop("Unexpected TranslationAI ORF format: missing tab-delimited coordinate payload")
    }

    transcript_id <- trimws(parts[[1]])
    coord_parts <- strsplit(trimws(parts[[2]]), ",")[[1]]
    coord_parts <- coord_parts[nzchar(coord_parts)]

    if (length(coord_parts) < 4) {
      stop("Unexpected TranslationAI ORF format: expected start, stop, and two scores")
    }

    data.frame(
      transcript_id = transcript_id,
      tis = as.numeric(coord_parts[[1]]),
      tts = as.numeric(coord_parts[[2]]),
      tis_score = as.numeric(coord_parts[[3]]),
      tts_score = as.numeric(coord_parts[[4]]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, parsed)
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

# Brief: Join TranslationAI predictions with transcript mapping for one selected gene.
# Args:
#   gene_query: Gene symbol or gene ID requested by the user.
#   orf_df: Parsed TranslationAI ORF prediction table.
#   mapping_df: Loaded transcript-to-gene mapping data frame.
# Returns: Prepared gene tis tts table object/table for the next step.
prepare_gene_tis_tts_table <- function(gene_query, orf_df, mapping_df) {
  gene_tx <- resolve_gene_transcripts(gene_query, mapping_df)
  merged_df <- gene_tx %>%
    inner_join(orf_df, by = "transcript_id") %>%
    arrange(tis, tts, transcript_id)

  if (nrow(merged_df) == 0) {
    stop(paste("No TranslationAI predictions matched the transcripts for gene:", gene_query))
  }

  merged_df$transcript_display <- factor(
    merged_df$transcript_id,
    levels = rev(unique(merged_df$transcript_id))
  )

  merged_df
}

# Brief: Write a gene-specific TIS/TTS coordinate table.
# Args:
#   gene_df: Gene-specific transcript table with prediction coordinates or scores.
#   output_file: Destination path for the generated table or figure.
# Returns: Writes gene-specific TIS/TTS CSV table.
write_gene_tis_tts_table <- function(gene_df, output_file) {
  ensure_dir(dirname(output_file))
  write.csv(gene_df, output_file, row.names = FALSE)
}

# Brief: Plot predicted TIS and TTS coordinates for transcripts of one gene.
# Args:
#   gene_df: Gene-specific transcript table with prediction coordinates or scores.
#   gene_query: Gene symbol or gene ID requested by the user.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for gene tis tts; writes a PDF when an output path is provided.
plot_gene_tis_tts <- function(gene_df,
                              gene_query,
                              output_file = NULL) {
  plot_df <- gene_df %>%
    mutate(
      transcript_display = factor(
        transcript_display,
        levels = levels(transcript_display)
      ),
      transcript_width = pmax(tts - tis, 1),
      tis_label = paste0("TIS ", tis, "\nscore=", sprintf("%.3f", tis_score)),
      tts_label = paste0("TTS ", tts, "\nscore=", sprintf("%.3f", tts_score))
    )

  x_min <- min(plot_df$tis, na.rm = TRUE)
  x_max <- max(plot_df$tts, na.rm = TRUE)
  x_span <- max(x_max - x_min, 1)
  x_pad <- max(60, 0.12 * x_span)

  p <- ggplot(plot_df, aes(y = transcript_display)) +
    geom_segment(
      aes(x = x_min - x_pad * 0.25, xend = x_max + x_pad * 0.25, yend = transcript_display),
      linewidth = 3.8,
      color = "#EEF3F8",
      lineend = "round"
    ) +
    geom_segment(
      aes(x = tis, xend = tts, yend = transcript_display),
      linewidth = 2.4,
      color = "#4C78A8",
      lineend = "round"
    ) +
    geom_point(
      aes(x = tis, fill = "TIS"),
      shape = 21,
      color = "white",
      stroke = 0.5,
      size = 3.4
    ) +
    geom_point(
      aes(x = tts, fill = "TTS"),
      shape = 21,
      color = "white",
      stroke = 0.5,
      size = 3.4
    ) +
    geom_text(
      aes(x = tis, label = tis_label),
      hjust = 1.02,
      nudge_x = -0.012 * x_span,
      vjust = 0.5,
      size = 3.1,
      lineheight = 0.92,
      color = "#C73E1D"
    ) +
    geom_text(
      aes(x = tts, label = tts_label),
      hjust = -0.02,
      nudge_x = 0.012 * x_span,
      vjust = 0.5,
      size = 3.1,
      lineheight = 0.92,
      color = "#1B8A5A"
    ) +
    scale_fill_manual(
      values = c(TIS = "#D95F02", TTS = "#1B9E77"),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(x_min - x_pad, x_max + x_pad),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Predicted TIS/TTS for ", gene_query),
      x = "Transcript coordinate",
      y = "Transcript"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2933"),
      axis.title = element_text(color = "#243447", face = "bold"),
      axis.text = element_text(color = "#243447"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#D9E2EC", linewidth = 0.35),
      legend.position = "top",
      legend.text = element_text(color = "#243447"),
      plot.margin = margin(12, 24, 12, 24)
    )

  if (!is.null(output_file) && output_file != "") {
    plot_height <- max(6, 0.85 * nrow(gene_df) + 2.8)
    # Keep gene panels readable but avoid overly wide PDFs.
    plot_width <- max(9.5, min(12, 9.5 + x_span / 4000))
    save_plot_pdf(p, output_file, width = plot_width, height = plot_height)
  }

  p
}
