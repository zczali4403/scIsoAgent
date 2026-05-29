# ------------------------------------------------------------------------------
# Script: exploring_isoforms_helpers.r
# Purpose: Helper utilities for exploratory isoform-level downstream analyses.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(data.table)
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
save_plot_pdf <- function(plot_obj, file_path, width = 8, height = 6) {
  ensure_dir(dirname(file_path))
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Read the transcript-to-gene mapping table produced during resource setup.
# Args:
#   mapping_file: Transcript-to-gene mapping CSV produced by resource setup.
# Returns: Loaded transcript gene mapping object/table from disk or package resources.
read_transcript_gene_mapping <- function(mapping_file) {
  mapping_df <- fread(mapping_file)

  required_cols <- c("transcript_id", "gene_id", "gene_symbol")
  missing_cols <- required_cols[!required_cols %in% colnames(mapping_df)]
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "The mapping file is missing required columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  mapping_df %>%
    dplyr::select(transcript_id, gene_id, gene_symbol) %>%
    dplyr::distinct(transcript_id, .keep_all = TRUE)
}

# Brief: Read a gene-id-to-gene-symbol mapping derived from the transcript mapping table.
# Args:
#   mapping_file: Transcript-to-gene mapping CSV produced by resource setup.
# Returns: Loaded gene symbol mapping object/table from disk or package resources.
read_gene_symbol_mapping <- function(mapping_file) {
  read_transcript_gene_mapping(mapping_file) %>%
    dplyr::select(gene_id, gene_symbol) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}

# Brief: Load a Seurat object and ensure a celltype metadata field is available for exploration.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   celltype_col: Metadata column containing cell type labels.
# Returns: Loaded Seurat object with a usable celltype metadata field for isoform exploration.
read_isoform_exploration_object <- function(input_rds, celltype_col = NULL) {
  seurat_obj <- readRDS(input_rds)
  if (!is.null(celltype_col) && celltype_col != "" && celltype_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- seurat_obj@meta.data[[celltype_col]]
  } else if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- Idents(seurat_obj)
  }
  seurat_obj
}

# Brief: Aggregate isoform expression by metadata group using Seurat pseudobulk aggregation.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
# Returns: Pseudobulk isoform expression data frame with feature IDs preserved as a column.
aggregate_isoform_expression <- function(seurat_obj,
                                         group_by = "celltype",
                                         assay = "isoform") {
  counts <- AggregateExpression(
    seurat_obj,
    assays = assay,
    return.seurat = FALSE,
    group.by = group_by
  )

  df <- as.data.frame(counts)
  df$feature_id <- rownames(df)
  df
}

# Brief: Join isoform features to transcript, gene ID, and gene-symbol annotations.
# Args:
#   pseudobulk_df: Pseudobulk expression data frame with features as rows.
#   transcript_gene_mapping: Mapping table with transcript_id, gene_id, and gene_symbol columns.
#   feature_col: Column containing feature identifiers in a table.
# Returns: Isoform expression table joined to transcript_id, gene_id, and gene_symbol annotations.
split_isoform_gene_ids <- function(pseudobulk_df,
                                   transcript_gene_mapping,
                                   feature_col = "feature_id") {
  if (missing(transcript_gene_mapping) || is.null(transcript_gene_mapping)) {
    stop("transcript_gene_mapping is required. Use read_transcript_gene_mapping() with the resource-setup dictionary.")
  }

  if (!feature_col %in% colnames(pseudobulk_df)) {
    stop(sprintf("feature column '%s' not found in pseudobulk_df", feature_col))
  }

  if (!all(c("transcript_id", "gene_id") %in% colnames(transcript_gene_mapping))) {
    stop("transcript_gene_mapping must include transcript_id and gene_id columns")
  }

  out <- pseudobulk_df %>%
    dplyr::rename(transcript_id = !!sym(feature_col)) %>%
    dplyr::left_join(
      transcript_gene_mapping %>% dplyr::select(transcript_id, gene_id, gene_symbol),
      by = "transcript_id"
    )

  missing_n <- sum(is.na(out$gene_id))
  if (missing_n > 0) {
    warning(sprintf("%d transcript features have no transcript→gene mapping and will be removed", missing_n))
    out <- out %>% dplyr::filter(!is.na(gene_id))
  }

  out
}

# Brief: Count distinct isoforms assigned to each gene.
# Args:
#   pseudobulk_df: Pseudobulk expression data frame with features as rows.
# Returns: Computed isoform count per gene table or vector.
compute_isoform_count_per_gene <- function(pseudobulk_df) {
  pseudobulk_df %>%
    group_by(gene_id) %>%
    summarise(n_isoforms = n_distinct(transcript_id), .groups = "drop")
}

# Brief: Assign genes to isoform-count bins such as 1, 2-3, 4-5, or 6+.
# Args:
#   isoform_count_per_gene: Table with one row per gene and its detected isoform count.
# Returns: Gene-level isoform-count table with an isoform_category column added.
categorize_isoform_count_per_gene <- function(isoform_count_per_gene) {
  isoform_count_per_gene %>%
    mutate(
      isoform_category = case_when(
        n_isoforms == 1 ~ "1",
        n_isoforms >= 2 & n_isoforms <= 3 ~ "2-3",
        n_isoforms >= 4 & n_isoforms <= 5 ~ "4-5",
        n_isoforms >= 6 ~ "6+"
      )
    )
}

# Brief: Summarize the number and percentage of genes in each isoform-count bin.
# Args:
#   categorized_df: Gene-level isoform-count table with assigned isoform-count categories.
# Returns: Summary table for isoform categories.
summarize_isoform_categories <- function(categorized_df) {
  categorized_df %>%
    dplyr::count(isoform_category) %>%
    mutate(percent = (n / sum(n)) * 100)
}

# Brief: Plot the distribution of isoform counts per gene.
# Args:
#   summary_df: Summary table used for plotting or writing results.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for isoform count summary; writes a PDF when an output path is provided.
plot_isoform_count_summary <- function(summary_df, output_file = NULL) {
  p <- ggplot(summary_df, aes(x = isoform_category, y = percent)) +
    geom_bar(stat = "identity", fill = "lightblue", color = "black") +
    labs(
      title = "Number of Isoforms per Gene",
      x = "Isoforms per Gene",
      y = "Genes, %"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 7, height = 5)
  }

  p
}

# Brief: Rank genes by the number of distinct isoforms detected.
# Args:
#   pseudobulk_df: Pseudobulk expression data frame with features as rows.
#   gene_symbol_mapping: Optional table mapping gene IDs to gene symbols for labels.
#   output_file: Destination path for the generated table or figure.
# Returns: Ranked genes by isoform count table.
rank_genes_by_isoform_count <- function(pseudobulk_df,
                                        gene_symbol_mapping = NULL,
                                        output_file = NULL) {
  gene_transcript_counts <- pseudobulk_df %>%
    group_by(gene_id) %>%
    summarise(unique_transcripts = n_distinct(transcript_id), .groups = "drop") %>%
    arrange(desc(unique_transcripts))

  if (!is.null(gene_symbol_mapping)) {
    gene_transcript_counts <- gene_transcript_counts %>%
      left_join(gene_symbol_mapping, by = "gene_id") %>%
      mutate(
        gene_symbol = ifelse(
          is.na(gene_symbol) | gene_symbol == "",
          gene_id,
          gene_symbol
        )
      )
  }

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    write.csv(gene_transcript_counts, output_file, row.names = FALSE)
  }

  gene_transcript_counts
}

# Brief: Plot the ranked isoform-rich gene distribution and label top genes.
# Args:
#   gene_transcript_counts: Table ranking genes by detected transcript count.
#   top_n: Number of top isoform-rich genes to label in the rank plot.
#   use_gene_symbol_labels: Whether to label highlighted genes by gene_symbol when available.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for top isoform rich genes; writes a PDF when an output path is provided.
plot_top_isoform_rich_genes <- function(gene_transcript_counts,
                                        top_n = 10,
                                        use_gene_symbol_labels = TRUE,
                                        output_file = NULL) {
  top_genes <- gene_transcript_counts %>% slice_max(unique_transcripts, n = top_n)

  plot_df <- gene_transcript_counts %>%
    arrange(desc(unique_transcripts)) %>%
    mutate(
      rank = row_number(),
      is_top = gene_id %in% top_genes$gene_id,
      label_text = if (use_gene_symbol_labels && "gene_symbol" %in% colnames(gene_transcript_counts)) {
        ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol)
      } else {
        gene_id
      }
    )

  p <- ggplot(plot_df, aes(x = rank, y = unique_transcripts)) +
    geom_point(color = "#6C8EBF", size = 0.9, alpha = 0.45) +
    geom_point(
      data = subset(plot_df, is_top),
      color = "#D95F02",
      size = 2.2,
      alpha = 0.95
    ) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Unique Transcripts per Gene",
      x = "Rank (log scale)",
      y = "# Transcripts (log scale)"
    ) +
    ggrepel::geom_label_repel(
      data = subset(plot_df, is_top),
      aes(label = label_text),
      fill = "#FFF7EC",
      color = "#2F2F2F",
      label.size = 0.25,
      size = 3,
      box.padding = 0.25,
      point.padding = 0.2,
      segment.color = "#8C8C8C",
      max.overlaps = Inf
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#222222"),
      axis.text = element_text(size = 10, color = "#333333"),
      axis.title = element_text(size = 11, color = "#222222"),
      panel.border = element_rect(color = "#444444", fill = NA, linewidth = 0.8),
      legend.position = "none"
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 8, height = 6)
  }

  p
}

# Brief: Rank one gene’s isoforms by total expression across cells.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   gene: Gene symbol or gene ID to analyze.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   slot: Assay data slot or layer to read from.
# Returns: Ranked isoforms for gene table.
rank_isoforms_for_gene <- function(seurat_obj,
                                   gene,
                                   assay = "isoform",
                                   slot = "data") {
  expression_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  matching_features <- grep(
    paste0("(^|-|\\b)", gene, "($|\\b)"),
    rownames(expression_matrix),
    value = TRUE
  )

  if (length(matching_features) == 0) {
    stop(paste("No isoform features matched the requested gene:", gene))
  }

  subset_expression <- expression_matrix[matching_features, , drop = FALSE]
  total_expression <- Matrix::rowSums(subset_expression)

  data.frame(
    feature = names(sort(total_expression, decreasing = TRUE)),
    total_expression = sort(total_expression, decreasing = TRUE),
    stringsAsFactors = FALSE
  )
}

# Brief: Write a selected gene’s ranked isoform table to disk.
# Args:
#   ranked_df: Ranked isoform or gene table to write or plot.
#   output_file: Destination path for the generated table or figure.
# Returns: Writes the selected gene isoform ranking CSV.
write_ranked_isoforms_for_gene <- function(ranked_df, output_file) {
  ensure_dir(dirname(output_file))
  write.csv(ranked_df, output_file, row.names = FALSE)
}

# Brief: Plot selected isoforms on a dimensional reduction.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
#   reduction: Dimensional reduction name used for embedding plots.
#   output_file: Destination path for the generated table or figure.
#   pt.size: Point size used in feature plots.
#   label.size: Label size used in feature plots.
# Returns: Plot object for selected isoforms featureplot; writes a PDF when an output path is provided.
plot_selected_isoforms_featureplot <- function(seurat_obj,
                                               features,
                                               reduction = "umap",
                                               output_file = NULL,
                                               pt.size = 0.4,
                                               label.size = 0.8) {
  plots <- lapply(features, function(one_feature) {
    p <- FeatureDimPlot(
      srt = seurat_obj,
      features = one_feature,
      reduction = reduction,
      theme_use = "theme_blank",
      pt.size = pt.size,
      label.size = label.size,
    )
    p <- p + theme(plot.title = element_text(size = 8))
    p
  })

  ncol_use <- min(3, length(plots))
  combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = ncol_use)

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(combined_plot, output_file, width = 16, height = 10)
  }

  combined_plot
}

# Brief: Plot selected isoform expression summaries across cell groups.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for selected isoforms stat; writes a PDF when an output path is provided.
plot_selected_isoforms_stat <- function(seurat_obj,
                                        features,
                                        group_by = "celltype",
                                        output_file = NULL) {
  plots <- lapply(features, function(one_feature) {
    p <- FeatureStatPlot(
      srt = seurat_obj,
      group.by = group_by,
      bg.by = group_by,
      stat.by = one_feature,
      add_box = TRUE
    )
    # Avoid duplicate-looking titles from multiple title layers.
    p <- p +
      labs(title = NULL, subtitle = NULL, caption = NULL) +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
    p
  })

  ncol_use <- min(3, length(plots))
  p <- cowplot::plot_grid(plotlist = plots, ncol = ncol_use)

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 16, height = 10)
  }

  p
}

# Brief: Generate a grouped heatmap for selected isoforms.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   features: Feature IDs, gene names, or isoform IDs to plot or summarize.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for selected isoforms heatmap; writes a PDF when an output path is provided.
plot_selected_isoforms_heatmap <- function(seurat_obj,
                                           features,
                                           group_by = "celltype",
                                           output_file = NULL) {
  ht <- GroupHeatmap(
    srt = seurat_obj,
    features = features,
    group.by = c(group_by),
    heatmap_palette = "YlOrRd",
    show_row_names = TRUE,
    row_names_side = "left",
    # Avoid duplicated transcript labels across composite heatmap panels.
    add_dot = FALSE,
    add_reticle = FALSE
  )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(ht$plot, output_file, width = 10, height = 8)
  }

  ht
}

# Brief: Write pseudobulk isoform expression in a format suitable for IsoVis input.
# Args:
#   pseudobulk_df: Pseudobulk expression data frame with features as rows.
#   output_file: Destination path for the generated table or figure.
# Returns: Writes the IsoVis pseudobulk CSV and invisibly returns the output path.
export_isovis_pseudobulk <- function(pseudobulk_df, output_file) {
  ensure_dir(dirname(output_file))
  export_df <- pseudobulk_df
  rownames(export_df) <- NULL
  write.csv(export_df, output_file, row.names = FALSE)
  invisible(output_file)
}
