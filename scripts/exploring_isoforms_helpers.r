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

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_plot_pdf <- function(plot_obj, file_path, width = 8, height = 6) {
  ensure_dir(dirname(file_path))
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

read_gene_symbol_mapping <- function(mapping_file) {
  mapping_df <- fread(mapping_file)

  required_cols <- c("gene_id", "gene_symbol")
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
    dplyr::select(gene_id, gene_symbol) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}

read_isoform_exploration_object <- function(input_rds, celltype_col = NULL) {
  seurat_obj <- readRDS(input_rds)
  if (!is.null(celltype_col) && celltype_col != "" && celltype_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- seurat_obj@meta.data[[celltype_col]]
  } else if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- Idents(seurat_obj)
  }
  seurat_obj
}

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

split_isoform_gene_ids <- function(pseudobulk_df,
                                   feature_col = "feature_id") {
  pseudobulk_df %>%
    tidyr::separate(
      !!sym(feature_col),
      into = c("transcript_id", "gene_id"),
      sep = "-(?=[^-]+$)"
    )
}

compute_isoform_count_per_gene <- function(pseudobulk_df) {
  pseudobulk_df %>%
    group_by(gene_id) %>%
    summarise(n_isoforms = n_distinct(transcript_id), .groups = "drop")
}

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

summarize_isoform_categories <- function(categorized_df) {
  categorized_df %>%
    dplyr::count(isoform_category) %>%
    mutate(percent = (n / sum(n)) * 100)
}

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

write_ranked_isoforms_for_gene <- function(ranked_df, output_file) {
  ensure_dir(dirname(output_file))
  write.csv(ranked_df, output_file, row.names = FALSE)
}

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
    save_plot_pdf(combined_plot, output_file, width = 16, height = 16)
  }

  combined_plot
}

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
        strip.text = element_text(size = 8)
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

export_isovis_pseudobulk <- function(pseudobulk_df, output_file) {
  ensure_dir(dirname(output_file))
  export_df <- pseudobulk_df
  rownames(export_df) <- NULL
  write.csv(export_df, output_file, row.names = FALSE)
  invisible(output_file)
}
