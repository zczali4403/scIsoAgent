# ------------------------------------------------------------------------------
# Script: isoform_classification_helpers.r
# Purpose: Helper utilities for SQANTI-based isoform classification summaries.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
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

# Brief: Select or recycle a project color palette for categorical plots.
# Args:
#   n: Number of items to select or display.
# Returns: Character vector of colors sampled or recycled from the project palette.
sample_palette_21 <- function(n) {
  if (n <= length(color_palette_21)) {
    sample(color_palette_21, n)
  } else {
    sample(rep(color_palette_21, length.out = n), n)
  }
}

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

# Brief: Read the SQANTI3 classification table.
# Args:
#   classification_file: Filesystem input or output location used by this step.
# Returns: Loaded sqanti classification object/table from disk or package resources.
read_sqanti_classification <- function(classification_file) {
  fread(classification_file, sep = "\t", data.table = FALSE)
}

# Brief: Load the Seurat object used for SQANTI summary analyses.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
#   group_col: Metadata column defining the groups for plotting or summaries.
# Returns: Loaded Seurat object for isoform classification summaries.
read_isoform_classification_object <- function(input_rds, group_col = NULL) {
  seurat_obj <- readRDS(input_rds)
  if (!is.null(group_col) && group_col != "" && group_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[[group_col]] <- seurat_obj@meta.data[[group_col]]
  }
  seurat_obj
}

# Brief: Documented helper for this workflow step.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
# Returns: Isoform pseudobulk expression data frame grouped by the requested metadata fields.
aggregate_isoform_pseudobulk <- function(seurat_obj,
                                         group_by = "sctype_db",
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

# Brief: Split isoform feature IDs into transcript ID and gene ID components.
# Args:
#   pseudobulk_df: Pseudobulk expression data frame with features as rows.
#   feature_col: Column containing feature identifiers in a table.
# Returns: Pseudobulk table with feature IDs separated into transcript_id and gene_id columns.
split_isoform_feature_ids <- function(pseudobulk_df, feature_col = "feature_id") {
  pseudobulk_df %>%
    tidyr::separate(
      !!sym(feature_col),
      into = c("transcript_id", "gene_id"),
      sep = "-(?=[^-]+$)"
    )
}

# Brief: Remove helper-generated prefixes from celltype labels.
# Args:
#   x: Celltype label string to clean.
# Returns: Normalized celltype label values.
normalize_celltype_label <- function(x) {
  # remove repeated helper prefixes like iso., iso.isoform., isoform.
  gsub("^(iso\\.|isoform\\.)+", "", x)
}

# Brief: Convert pseudobulk isoform expression into long format by cell type.
# Args:
#   pseudobulk_data: Pseudobulk isoform expression table with feature annotations.
#   group_prefix: Column-name prefix used to identify grouped pseudobulk expression columns.
# Returns: Constructed long isoform expression object/table.
build_long_isoform_expression <- function(pseudobulk_data, group_prefix = "iso.") {
  pseudobulk_data %>%
    pivot_longer(
      cols = starts_with(group_prefix),
      names_to = "cell_type",
      values_to = "expression"
    ) %>%
    mutate(cell_type = normalize_celltype_label(cell_type))
}

# Brief: Count expressed genes and isoforms per cell type.
# Args:
#   long_data: Long-format expression table used for grouped summaries.
# Returns: Summary table for celltype gene isoform counts.
summarize_celltype_gene_isoform_counts <- function(long_data) {
  filtered_data <- long_data %>%
    filter(expression > 0) %>%
    dplyr::select(cell_type, transcript_id, gene_id) %>%
    distinct()

  filtered_data %>%
    group_by(cell_type) %>%
    summarise(
      num_genes = n_distinct(gene_id),
      num_isoforms = n_distinct(transcript_id),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(num_genes, num_isoforms),
      names_to = "Category",
      values_to = "Count"
    )
}

# Brief: Plot gene and isoform counts across cell types.
# Args:
#   cell_type_summary: Celltype-level gene/isoform count summary table.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for celltype gene isoform counts; writes a PDF when an output path is provided.
plot_celltype_gene_isoform_counts <- function(cell_type_summary, output_file = NULL) {
  fill_levels <- unique(as.character(cell_type_summary$Category))
  fill_values <- stats::setNames(sample_palette_21(length(fill_levels)), fill_levels)

  p <- ggplot(cell_type_summary, aes(x = cell_type, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.8) +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    labs(
      title = "Number of genes and isoforms \n per cell type",
      x = "Cell Type",
      y = "Count"
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_blank()
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 9, height = 6)
  }

  p
}

# Brief: Join SQANTI classifications onto pseudobulk isoform expression records.
# Args:
#   pseudobulk_data: Pseudobulk isoform expression table with feature annotations.
#   sqanti_df: SQANTI3 classification table to join onto isoform expression data.
# Returns: Merged sqanti with pseudobulk object/table.
merge_sqanti_with_pseudobulk <- function(pseudobulk_data, sqanti_df) {
  merge(
    pseudobulk_data,
    sqanti_df,
    by.x = "transcript_id",
    by.y = "isoform",
    all.x = TRUE
  )
}

# Brief: Summarize SQANTI structural categories, subcategories, and coding status by cell type.
# Args:
#   merged_data: Pseudobulk isoform expression table already merged with SQANTI annotations.
#   group_prefix: Column-name prefix used to identify grouped pseudobulk expression columns.
# Returns: List with long expression data and structural/coding category summary table.
summarize_sqanti_categories_by_celltype <- function(merged_data, group_prefix = "iso.") {
  long_data <- merged_data %>%
    pivot_longer(
      cols = starts_with(group_prefix),
      names_to = "cell_type",
      values_to = "expression"
    ) %>%
    mutate(cell_type = normalize_celltype_label(cell_type)) %>%
    filter(expression > 0) %>%
    distinct()

  category_summary <- long_data %>%
    group_by(cell_type, structural_category, subcategory, coding) %>%
    summarise(
      num_genes = n_distinct(gene_id),
      num_isoforms = n_distinct(transcript_id),
      .groups = "drop"
    )

  list(
    long_data = long_data,
    category_summary = category_summary
  )
}

# Brief: Plot isoform structural-category counts across cell types.
# Args:
#   category_summary: SQANTI category summary table by cell type.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for structural category summary; writes a PDF when an output path is provided.
plot_structural_category_summary <- function(category_summary, output_file = NULL) {
  fill_levels <- unique(as.character(category_summary$structural_category))
  fill_values <- stats::setNames(sample_palette_21(length(fill_levels)), fill_levels)

  p <- ggplot(category_summary, aes(x = cell_type, y = num_isoforms, fill = structural_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    labs(
      title = "Isoforms by structural category \n across cell types",
      x = "Cell Type",
      y = "Number of Isoforms",
      fill = "Structural Category"
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 10, height = 7)
  }

  p
}

# Brief: Plot SQANTI structural-subcategory counts across cell types.
# Args:
#   category_summary: SQANTI category summary table by cell type.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for structural subcategory summary; writes a PDF when an output path is provided.
plot_structural_subcategory_summary <- function(category_summary, output_file = NULL) {
  fill_levels <- unique(as.character(category_summary$subcategory))
  fill_values <- stats::setNames(sample_palette_21(length(fill_levels)), fill_levels)

  p <- ggplot(category_summary, aes(x = cell_type, y = num_isoforms, fill = subcategory)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    labs(
      title = "Isoforms by structural subcategory \n across cell types",
      x = "Cell Type",
      y = "Number of Isoforms",
      fill = "subcategory"
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 10, height = 7)
  }

  p
}

# Brief: Plot the coding versus noncoding isoform proportions across cell types.
# Args:
#   category_summary: SQANTI category summary table by cell type.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for coding proportion summary; writes a PDF when an output path is provided.
plot_coding_proportion_summary <- function(category_summary, output_file = NULL) {
  fill_levels <- unique(as.character(category_summary$coding))
  fill_values <- stats::setNames(sample_palette_21(length(fill_levels)), fill_levels)

  p <- ggplot(category_summary, aes(x = cell_type, y = num_isoforms, fill = coding)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    labs(
      title = "Proportion of coding isoforms \n across cell types",
      x = "Cell Type",
      y = "Proportion of Isoforms",
      fill = "Structural Category"
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 10, height = 7)
  }

  p
}

# Brief: Combine major SQANTI summary plots into one overview figure.
# Args:
#   cell_type_summary: Celltype-level gene/isoform count summary table.
#   category_summary: SQANTI category summary table by cell type.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for sqanti classification overview; writes a PDF when an output path is provided.
plot_sqanti_classification_overview <- function(cell_type_summary,
                                                category_summary,
                                                output_file = NULL) {
  p1 <- plot_celltype_gene_isoform_counts(cell_type_summary)
  p2 <- plot_structural_category_summary(category_summary)
  p3 <- plot_structural_subcategory_summary(category_summary)
  p4 <- plot_coding_proportion_summary(category_summary)

  p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 16, height = 12)
  }

  p
}

# Brief: Identify isoforms expressed above cutoff in only one cell type.
# Args:
#   long_data: Long-format expression table used for grouped summaries.
#   expression_cutoff: Expression cutoff used to define celltype-specific isoform presence.
#   max_expression_in_all_other_cells_types: Maximum allowed expression in non-target cell types for specificity calls.
# Returns: Long-format rows for isoforms passing the selected celltype-specific expression rule.
find_celltype_specific_isoforms <- function(long_data,
                                            expression_cutoff = 5,
                                            max_expression_in_all_other_cells_types = 1) {
  long_data %>%
    group_by(transcript_id) %>%
    filter(
      sum(expression > expression_cutoff) <= 1,
      all(expression < max_expression_in_all_other_cells_types | expression > expression_cutoff)
    ) %>%
    ungroup()
}

# Brief: Count celltype-specific isoforms for each cell type.
# Args:
#   exclusive_isoforms: Transcript/isoform identifier or transcript-level data used by the analysis.
#   expression_cutoff: Expression cutoff used to define celltype-specific isoform presence.
# Returns: Table of unique celltype-specific isoform counts by cell type.
count_celltype_specific_isoforms <- function(exclusive_isoforms,
                                             expression_cutoff = 5) {
  exclusive_isoforms %>%
    group_by(cell_type) %>%
    summarise(
      unique_isoforms = sum(expression > expression_cutoff),
      .groups = "drop"
    )
}

# Brief: Plot the number of celltype-specific isoforms per cell type.
# Args:
#   exclusive_isoforms_count: Transcript/isoform identifier or transcript-level data used by the analysis.
#   output_file: Destination path for the generated table or figure.
# Returns: Plot object for celltype specific isoform counts; writes a PDF when an output path is provided.
plot_celltype_specific_isoform_counts <- function(exclusive_isoforms_count,
                                                  output_file = NULL) {
  fill_levels <- unique(as.character(exclusive_isoforms_count$cell_type))
  fill_values <- stats::setNames(sample_palette_21(length(fill_levels)), fill_levels)

  p <- ggplot(exclusive_isoforms_count, aes(x = cell_type, y = unique_isoforms, fill = cell_type)) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    labs(
      title = "Number of isoforms unique \n to each cell type",
      x = "Cell Type",
      y = "Number of Unique Isoforms"
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_blank()
    )

  if (!is.null(output_file) && output_file != "") {
    save_plot_pdf(p, output_file, width = 9, height = 6)
  }

  p
}
