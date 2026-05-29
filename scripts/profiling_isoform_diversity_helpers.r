# ------------------------------------------------------------------------------
# Script: profiling_isoform_diversity_helpers.r
# Purpose: Helper utilities for entropy-based isoform diversity profiling and plotting.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(data.table)
  library(pheatmap)
  library(ggbeeswarm)
  library(cowplot)
  library(scales)
})

color_palette_21 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#00CED1", "#A65628", "#F781BF", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#8B4513", "#FFFF33", "#FF69B4", "#2E8B57",
  "#C8A2C8"
)

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

# Brief: Locate the isoform-to-gene-symbol dictionary for diversity labeling.
# Args:
#   gene_symbol_dict_csv: Optional path to the resource-setup isoform_gene_symbol_dict.csv file.
#   output_rds: Path where the output Seurat RDS object should be saved.
#   output_root: Workflow results root used to auto-discover the resource-setup dictionary.
# Returns: Resolved gene symbol dict value or matching records.
resolve_gene_symbol_dict <- function(gene_symbol_dict_csv = NULL, output_rds = NULL, output_root = NULL) {
  if (!is.null(gene_symbol_dict_csv) && gene_symbol_dict_csv != "" && file.exists(gene_symbol_dict_csv)) {
    return(gene_symbol_dict_csv)
  }

  if (is.null(output_root) || output_root == "") {
    if (!is.null(output_rds) && output_rds != "") {
      output_root <- dirname(dirname(dirname(output_rds)))
    }
  }

  if (!is.null(output_root) && output_root != "") {
    rs_dir <- file.path(output_root, "resource_setup")
    if (dir.exists(rs_dir)) {
      cands <- list.files(
        rs_dir,
        pattern = "isoform_gene_symbol_dict\\.csv$",
        full.names = TRUE
      )
      if (length(cands) > 0) {
        cands <- cands[order(nchar(cands))]
        return(cands[[1]])
      }
    }
  }

  NULL
}

# Brief: Load and validate a gene-symbol dictionary table.
# Args:
#   dict_csv: Candidate gene-symbol dictionary CSV path.
# Returns: Loaded gene symbol dict object/table from disk or package resources.
load_gene_symbol_dict <- function(dict_csv) {
  if (is.null(dict_csv) || dict_csv == "" || !file.exists(dict_csv)) {
    return(NULL)
  }

  dict <- data.table::fread(dict_csv, data.table = FALSE)
  req_cols <- c("gene_id", "gene_symbol")
  if (!all(req_cols %in% colnames(dict))) {
    warning("Gene symbol dictionary is missing required columns: gene_id/gene_symbol")
    return(NULL)
  }

  dict <- dict %>%
    dplyr::select(gene_id, gene_symbol) %>%
    dplyr::mutate(
      gene_id = as.character(gene_id),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    dplyr::filter(!is.na(gene_id), gene_id != "") %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)

  dict
}

# Brief: Replace gene IDs with gene symbols while preserving unique row names.
# Args:
#   gene_ids: Gene identifiers to map or summarize.
#   dict_df: Loaded gene ID to gene symbol dictionary table.
# Returns: Character vector of gene symbols with unmapped IDs preserved and names made unique.
map_gene_ids_to_symbols <- function(gene_ids, dict_df = NULL) {
  if (is.null(dict_df) || length(gene_ids) == 0) {
    return(gene_ids)
  }

  idx <- match(gene_ids, dict_df$gene_id)
  mapped <- dict_df$gene_symbol[idx]
  mapped[is.na(mapped) | mapped == ""] <- gene_ids[is.na(mapped) | mapped == ""]
  make.unique(mapped)
}

# Brief: Load the Seurat object used for entropy-based isoform diversity analysis.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   celltype_col: Metadata column containing cell type labels.
# Returns: Loaded Seurat object with a usable celltype metadata field for diversity analysis.
read_isoform_diversity_object <- function(input_rds, celltype_col = NULL) {
  seurat_obj <- readRDS(input_rds)
  if (!is.null(celltype_col) && celltype_col != "" && celltype_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- seurat_obj@meta.data[[celltype_col]]
  } else if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- Idents(seurat_obj)
  }
  seurat_obj
}

# Brief: Subset the Seurat object to one selected cell type for diversity analysis.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   subset_col: Metadata column used to choose a cell subset.
#   subset_value: Metadata value identifying the cell subset to keep.
# Returns: Seurat object containing only cells matching the requested metadata value.
subset_diversity_celltype <- function(seurat_obj, subset_col, subset_value) {
  if (!subset_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Metadata column not found: ", subset_col))
  }

  keep_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[subset_col]] == subset_value]
  if (length(keep_cells) == 0) {
    stop(paste0("No cells found for ", subset_col, " == ", subset_value))
  }

  subset(seurat_obj, cells = keep_cells)
}

# Brief: Convert Seurat isoform counts to SingleCellExperiment format and compute isoform entropy.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   genes: Gene identifiers or symbols selected for analysis.
#   min_counts_per_cell: Minimum gene or isoform counts per cell required for entropy analysis.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   output_rds: Path where the output Seurat RDS object should be saved.
#   threads: Number of worker threads to use.
#   show_progress: Whether to display progress messages during long computations.
#   gene_symbol_dict_csv: Optional path to the resource-setup isoform_gene_symbol_dict.csv file.
#   map_to_gene_symbol: Whether entropy result rows should be relabeled from gene IDs to gene symbols.
# Returns: List containing normalized_entropy and raw_entropy matrices; optionally saved as RDS.
run_isoform_diversity <- function(seurat_obj,
                                  genes = rownames(seurat_obj),
                                  min_counts_per_cell = 10,
                                  assay = "iso",
                                  output_rds = NULL,
                                  threads = 1,
                                  show_progress = interactive(),
                                  gene_symbol_dict_csv = NULL,
                                  map_to_gene_symbol = TRUE) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste0("Assay not found in Seurat object: ", assay))
  }

  assay_obj <- seurat_obj@assays[[assay]]

  # S4-style extraction of counts matrix (Assay5 or legacy Assay)
  if ("layers" %in% slotNames(assay_obj)) {
    if (!"counts" %in% names(assay_obj@layers)) {
      stop(paste0("Assay '", assay, "' does not contain a 'counts' layer"))
    }
    expr <- assay_obj@layers[["counts"]]
  } else if ("counts" %in% slotNames(assay_obj)) {
    expr <- assay_obj@counts
  } else {
    stop(paste0("Unable to locate counts matrix in assay: ", assay))
  }

  # Ensure sparse matrix and valid dimnames for downstream split()
  if (!inherits(expr, "sparseMatrix")) {
    expr <- as(expr, "CsparseMatrix")
  }

  feature_ids <- rownames(expr)
  if (is.null(feature_ids) || length(feature_ids) != nrow(expr)) {
    if ("features" %in% slotNames(assay_obj) && !is.null(rownames(assay_obj@features))) {
      feature_ids <- rownames(assay_obj@features)
      rownames(expr) <- feature_ids
    } else {
      stop("Counts matrix has no valid feature names (rownames)")
    }
  }

  cell_ids <- colnames(expr)
  if (is.null(cell_ids) || length(cell_ids) != ncol(expr)) {
    cell_ids <- colnames(seurat_obj)
    colnames(expr) <- cell_ids
  }

  if (!is.null(genes) && length(genes) > 0) {
    genes <- intersect(genes, rownames(expr))
    if (length(genes) == 0) {
      warning("Provided genes do not overlap with assay features; using all assay features instead")
      genes <- rownames(expr)
    }
    expr <- expr[genes, , drop = FALSE]
  }

  feature_ids <- rownames(expr)

  # Prefer transcript->gene mapping from resource-setup dictionary.
  resolved_dict <- resolve_gene_symbol_dict(
    gene_symbol_dict_csv = gene_symbol_dict_csv,
    output_rds = output_rds
  )
  gene_ids <- NULL
  if (!is.null(resolved_dict) && file.exists(resolved_dict)) {
    dict_full <- data.table::fread(resolved_dict, data.table = FALSE)
    if (all(c("transcript_id", "gene_id") %in% colnames(dict_full))) {
      dict_full <- dict_full %>%
        dplyr::select(transcript_id, gene_id) %>%
        dplyr::mutate(
          transcript_id = as.character(transcript_id),
          gene_id = as.character(gene_id)
        ) %>%
        dplyr::filter(!is.na(transcript_id), transcript_id != "", !is.na(gene_id), gene_id != "") %>%
        dplyr::distinct(transcript_id, .keep_all = TRUE)

      idx <- match(feature_ids, dict_full$transcript_id)
      gene_ids <- dict_full$gene_id[idx]

      # Fallback by stripping version suffix if needed (e.g., ENST... .1)
      miss <- is.na(gene_ids)
      if (any(miss)) {
        idx2 <- match(sub("\\.[0-9]+$", "", feature_ids[miss]), sub("\\.[0-9]+$", "", dict_full$transcript_id))
        gene_ids[miss] <- dict_full$gene_id[idx2]
      }
    }
  }

  # Final fallback: parse transcript-gene pattern from feature ids
  if (is.null(gene_ids)) {
    split_ids <- strsplit(feature_ids, "-(?=[^-]+$)", perl = TRUE)
    gene_ids <- as.character(vapply(split_ids, function(x) {
      if (length(x) >= 2) x[[2]] else x[[1]]
    }, character(1)))
  } else {
    # Keep unmapped features as-is to avoid missing rowData
    gene_ids[is.na(gene_ids) | gene_ids == ""] <- feature_ids[is.na(gene_ids) | gene_ids == ""]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = setNames(list(expr), assay)
  )
  SummarizedExperiment::rowData(sce)$gene_id <- gene_ids

  entropy_mat <- find_diversity(
    sce = sce,
    assay = assay,
    gene_col = "gene_id",
    min_counts_per_cell = min_counts_per_cell,
    threads = threads,
    show_progress = show_progress
  )

  if (map_to_gene_symbol) {
    resolved_dict <- resolve_gene_symbol_dict(
      gene_symbol_dict_csv = gene_symbol_dict_csv,
      output_rds = output_rds
    )
    dict_df <- load_gene_symbol_dict(resolved_dict)
    if (!is.null(dict_df)) {
      mapped_names <- map_gene_ids_to_symbols(rownames(entropy_mat), dict_df)
      rownames(entropy_mat) <- mapped_names
    } else {
      warning("No gene symbol dictionary resolved; keeping gene_id rownames")
    }
  }

  res <- list(
    normalized_entropy = entropy_mat,
    raw_entropy = entropy_mat
  )

  if (!is.null(output_rds) && output_rds != "") {
    ensure_dir(dirname(output_rds))
    saveRDS(res, output_rds)
  }

  res
}

# Brief: Draw a heatmap of normalized isoform entropy values.
# Args:
#   normalized_entropy: Matrix of normalized isoform entropy values.
#   output_file: Destination path for the generated table or figure.
#   cluster_rows: Whether heatmap rows should be clustered.
#   cluster_cols: Whether heatmap columns should be clustered.
#   main: Main title text for a heatmap.
#   na_col: Color used for missing entropy values in the heatmap.
# Returns: Plot object for entropy heatmap; writes a PDF when an output path is provided.
plot_entropy_heatmap <- function(normalized_entropy,
                                 output_file = NULL,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 main = "Entropy heatmap",
                                 na_col = "darkgrey") {
  palette_fun <- colorRampPalette(c("#3288BD", "white", "#D53E4F"))
  heat_colors <- palette_fun(100)

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    pdf(output_file, width = 10, height = 8, bg = "white")
    pheatmap::pheatmap(
      normalized_entropy,
      show_rownames = FALSE,
      show_colnames = FALSE,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      main = main,
      na_col = na_col,
      color = heat_colors
    )
    dev.off()
    return(invisible(NULL))
  }

  pheatmap::pheatmap(
    normalized_entropy,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    main = main,
    na_col = na_col,
    color = heat_colors
  )
}

# Brief: Compute median isoform entropy per gene.
# Args:
#   normalized_entropy: Matrix of normalized isoform entropy values.
# Returns: Summary table for entropy by gene.
summarize_entropy_by_gene <- function(normalized_entropy) {
  data.frame(
    gene = rownames(normalized_entropy),
    median_entropy = apply(normalized_entropy, 1, median, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Brief: Select genes with the lowest and highest median isoform entropy.
# Args:
#   entropy_summary: Per-gene entropy summary table.
#   n: Number of items to select or display.
# Returns: List containing global median plus low- and high-entropy gene tables.
get_low_high_entropy_genes <- function(entropy_summary, n = 25) {
  global_median <- median(entropy_summary$median_entropy, na.rm = TRUE)

  df_nonzero <- entropy_summary %>%
    dplyr::filter(median_entropy > 0) %>%
    dplyr::mutate(median_entropy = as.numeric(median_entropy))

  df_low <- df_nonzero %>%
    dplyr::arrange(median_entropy) %>%
    dplyr::slice_head(n = n)

  df_high <- df_nonzero %>%
    dplyr::arrange(desc(median_entropy)) %>%
    dplyr::slice_head(n = n)

  list(
    global_median = global_median,
    low = df_low,
    high = df_high
  )
}

# Brief: Plot high-entropy and low-entropy genes with a global median reference.
# Args:
#   entropy_summary: Per-gene entropy summary table.
#   n: Number of items to select or display.
#   output_file: Destination path for the generated table or figure.
# Returns: List containing combined plot, low/high gene tables, and global median.
plot_high_low_entropy_genes <- function(entropy_summary,
                                        n = 25,
                                        output_file = NULL) {
  split_res <- get_low_high_entropy_genes(entropy_summary, n = n)

  p_low <- ggplot(split_res$low, aes(x = median_entropy, y = reorder(gene, median_entropy), fill = median_entropy)) +
    geom_col() +
    geom_vline(xintercept = split_res$global_median, linetype = "dashed", color = "black") +
    scale_fill_gradient(low = "darkblue", high = "skyblue") +
    labs(x = "Median Normalised Entropy", y = "Low Entropy Genes") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  p_high <- ggplot(split_res$high, aes(x = median_entropy, y = reorder(gene, median_entropy), fill = median_entropy)) +
    geom_col() +
    geom_vline(xintercept = split_res$global_median, linetype = "dashed", color = "black") +
    scale_fill_gradient(low = "mistyrose", high = "firebrick") +
    labs(x = "Median Normalised Entropy", y = "High Entropy Genes") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  combined_plot <- cowplot::plot_grid(p_high, p_low, ncol = 2)

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(combined_plot, output_file, width = 14, height = 8)
  }

  list(
    plot = combined_plot,
    low = split_res$low,
    high = split_res$high,
    global_median = split_res$global_median
  )
}

# Brief: Reshape selected gene entropy values into long format with cell metadata.
# Args:
#   normalized_entropy: Matrix of normalized isoform entropy values.
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   genes: Gene identifiers or symbols selected for analysis.
#   celltype_col: Metadata column containing cell type labels.
# Returns: Prepared gene entropy long object/table for the next step.
prepare_gene_entropy_long <- function(normalized_entropy,
                                      seurat_obj,
                                      genes,
                                      celltype_col = "celltype") {
  df_entropy <- normalized_entropy[genes, , drop = FALSE]
  df_long <- df_entropy %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "cell",
      values_to = "entropy"
    )

  df_long$celltype <- seurat_obj@meta.data[df_long$cell, celltype_col]
  df_long
}

# Brief: Plot per-cell entropy values for one selected gene.
# Args:
#   df_long: Long-format entropy or expression table used for plotting.
#   gene_name: Gene label used in a gene-specific plot.
#   group_col: Metadata column defining the groups for plotting or summaries.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for gene entropy; writes a PDF when an output path is provided.
plot_gene_entropy <- function(df_long,
                              gene_name,
                              group_col = "celltype",
                              output_file = NULL) {
  df_gene <- df_long %>% dplyr::filter(gene == gene_name)
  med_val <- median(df_gene$entropy, na.rm = TRUE)

  p <- ggplot(df_gene, aes_string(x = group_col, y = "entropy")) +
    ggbeeswarm::geom_quasirandom(alpha = 0.4, varwidth = TRUE) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.6,
      fatten = 2,
      color = "red"
    ) +
    labs(
      title = paste0("Isoform entropy: ", gene_name),
      subtitle = paste0("Median = ", round(med_val, 3)),
      x = "Cell type",
      y = "Normalised Shannon entropy"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.subtitle = element_text(color = "grey40")
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 8, height = 6)
  }

  p
}

# Brief: Combine gene-specific entropy plots for multiple selected genes.
# Args:
#   df_long: Long-format entropy or expression table used for plotting.
#   genes: Gene identifiers or symbols selected for analysis.
#   group_col: Metadata column defining the groups for plotting or summaries.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for gene entropy batch; writes a PDF when an output path is provided.
plot_gene_entropy_batch <- function(df_long,
                                    genes,
                                    group_col = "celltype",
                                    output_file = NULL) {
  plots <- lapply(genes, function(gene_name) {
    plot_gene_entropy(df_long, gene_name, group_col = group_col)
  })
  combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2)

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(combined_plot, output_file, width = 14, height = 8)
  }

  combined_plot
}

# Brief: Draw per-cell pie charts of isoform proportions for one gene.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   gene: Gene symbol or gene ID to analyze.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   slot: Assay data slot or layer to read from.
#   top_n: Number of dominant isoforms retained separately before grouping the rest as other.
#   n_cells: Number of cells to include in a per-cell visualization.
#   plot_col: Number of columns in a multi-panel plot layout.
#   output_file: Destination path for the generated table or figure.
#   gene_symbol_dict_csv: Optional path to the resource-setup isoform_gene_symbol_dict.csv file.
#   output_root: Workflow results root used to auto-discover the resource-setup dictionary.
# Returns: Plot object for plotisoformpiecells; writes a PDF when an output path is provided.
plotIsoformPieCells <- function(seurat_obj,
                                gene,
                                assay = "iso",
                                slot = "counts",
                                top_n = 3,
                                n_cells = 12,
                                plot_col = 4,
                                output_file = NULL,
                                gene_symbol_dict_csv = NULL,
                                output_root = NULL) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste0("Assay not found in Seurat object: ", assay))
  }

  assay_obj <- seurat_obj@assays[[assay]]
  if ("layers" %in% slotNames(assay_obj)) {
    if (!slot %in% names(assay_obj@layers)) {
      stop(paste0("Layer not found in assay '", assay, "': ", slot))
    }
    expr <- assay_obj@layers[[slot]]
  } else if (slot %in% slotNames(assay_obj)) {
    expr <- slot(assay_obj, slot)
  } else {
    stop(paste0("Slot not found in assay '", assay, "': ", slot))
  }

  if (!inherits(expr, "sparseMatrix")) {
    expr <- as(expr, "CsparseMatrix")
  }

  feature_ids <- rownames(expr)
  if (is.null(feature_ids) || length(feature_ids) == 0) {
    stop("Assay matrix has no feature names")
  }

  # Parse feature identifiers
  transcript_ids <- sub("-([^-]+)$", "", feature_ids, perl = TRUE)
  parsed_gene_ids <- sub("^.*-([^-]+)$", "\\1", feature_ids, perl = TRUE)

  # Resolve dictionary for gene symbol -> gene/transcript mapping
  if ((is.null(output_root) || output_root == "") && !is.null(output_file) && output_file != "") {
    output_root <- dirname(dirname(dirname(output_file)))
  }

  dict_path <- resolve_gene_symbol_dict(
    gene_symbol_dict_csv = gene_symbol_dict_csv,
    output_root = output_root
  )

  target_gene_ids <- character(0)
  target_transcript_ids <- character(0)
  if (!is.null(dict_path) && file.exists(dict_path)) {
    dict_full <- data.table::fread(dict_path, data.table = FALSE)
    if ("gene_symbol" %in% colnames(dict_full)) {
      dict_full$gene_symbol <- as.character(dict_full$gene_symbol)
      sel <- dict_full$gene_symbol == gene
      if ("gene_id" %in% colnames(dict_full)) {
        target_gene_ids <- unique(as.character(dict_full$gene_id[sel]))
      }
      if ("transcript_id" %in% colnames(dict_full)) {
        target_transcript_ids <- unique(as.character(dict_full$transcript_id[sel]))
      }
    }
  }

  match_idx <- rep(FALSE, length(feature_ids))
  # Direct feature contains requested token (for gene_id input or literal match)
  match_idx <- match_idx | grepl(paste0("(^|[-_])", gene, "($|[-_])"), feature_ids)
  # Parsed gene_id/transcript_id matches
  match_idx <- match_idx | (parsed_gene_ids %in% c(gene, target_gene_ids))
  match_idx <- match_idx | (transcript_ids %in% target_transcript_ids)

  matching_features <- feature_ids[match_idx]
  if (length(matching_features) == 0) {
    stop(paste("No isoform features matched the requested gene:", gene))
  }

  gene_expr <- expr[matching_features, , drop = FALSE]
  cell_totals <- Matrix::colSums(gene_expr)
  n_keep <- min(n_cells, sum(cell_totals > 0))
  if (n_keep == 0) {
    stop("No cells with positive isoform counts were found for the requested gene")
  }
  selected_cells <- names(sort(cell_totals, decreasing = TRUE))[seq_len(n_keep)]

  top_features_global <- names(sort(Matrix::rowSums(gene_expr), decreasing = TRUE))[seq_len(min(top_n, nrow(gene_expr)))]
  fill_levels <- c(top_features_global, "other")
  palette_vals <- rep(color_palette_21, length.out = length(fill_levels))
  fill_values <- setNames(palette_vals, fill_levels)
  fill_values["other"] <- "#B3B3B3"

  plot_list <- lapply(selected_cells, function(cell_id) {
    cell_counts <- as.numeric(gene_expr[, cell_id])
    names(cell_counts) <- rownames(gene_expr)
    cell_counts <- sort(cell_counts, decreasing = TRUE)
    top_features_cell <- names(cell_counts)[seq_len(min(top_n, length(cell_counts)))]

    top_values <- cell_counts[top_features_cell]
    other_value <- sum(cell_counts[setdiff(names(cell_counts), top_features_cell)])

    pie_df <- data.frame(
      transcript = c(top_features_cell, "other"),
      count = c(top_values, other_value),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(count > 0)

    ggplot(pie_df, aes(x = "", y = count, fill = transcript)) +
      geom_col(width = 1, color = "white") +
      coord_polar(theta = "y") +
      scale_fill_manual(values = fill_values[names(fill_values) %in% pie_df$transcript]) +
      labs(title = cell_id, fill = NULL) +
      theme_void() +
      theme(
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.position = "none"
      )
  })

  pie_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = plot_col)

  legend_df <- data.frame(
    transcript = factor(fill_levels, levels = fill_levels),
    count = 1,
    stringsAsFactors = FALSE
  )

  legend_plot <- ggplot(legend_df, aes(x = "", y = count, fill = transcript)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = fill_values) +
    theme_void() +
    theme(legend.position = "right")

  shared_legend <- cowplot::get_legend(legend_plot)
  combined_plot <- cowplot::plot_grid(pie_grid, shared_legend, ncol = 2, rel_widths = c(5, 1))

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(combined_plot, output_file, width = 4 * plot_col + 2, height = 4 * ceiling(length(plot_list) / plot_col))
  }

  combined_plot
}
