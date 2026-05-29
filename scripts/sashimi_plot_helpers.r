# ------------------------------------------------------------------------------
# Script: sashimi_plot_helpers.r
# Purpose: Helper utilities for barcode extraction, BAM subsetting, and sashimi plotting.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(rtracklayer)
})

# Brief: Create an output directory when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the target directory exists; no analysis object is returned.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Brief: Load the annotated Seurat object used to extract barcode groups for sashimi plots.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
# Returns: Loaded Seurat object for sashimi barcode extraction.
read_sashimi_seurat_object <- function(input_rds) {
  readRDS(input_rds)
}

# Brief: List metadata columns available for grouping sashimi comparisons.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
# Returns: Character vector listing available metadata columns.
list_metadata_columns <- function(seurat_obj) {
  colnames(seurat_obj@meta.data)
}

# Brief: Summarize value counts for one metadata column.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   column: Metadata column whose values should be counted and reported.
#   max_values: Maximum number of most frequent metadata values to report.
# Returns: Summary table for metadata column.
summarize_metadata_column <- function(seurat_obj,
                                      column,
                                      max_values = 20) {
  if (!column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Metadata column not found:", column))
  }

  values <- as.character(seurat_obj@meta.data[[column]])
  summary_df <- data.frame(value = values, stringsAsFactors = FALSE) %>%
    filter(!is.na(value) & value != "") %>%
    count(value, name = "n", sort = TRUE)

  if (nrow(summary_df) > max_values) {
    summary_df <- summary_df[seq_len(max_values), , drop = FALSE]
  }

  summary_df
}

# Brief: Remove sample prefixes from Seurat barcodes before BAM tag matching.
# Args:
#   barcodes: Cell barcode or cell-level metadata used for filtering or grouping.
# Returns: Barcode vector with the leading sample prefix removed.
strip_barcode_prefix <- function(barcodes) {
  sub("^[^_]+_", "", barcodes)
}

# Brief: Build a barcode metadata table with sample, group, and optional celltype labels.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   sample_col: Metadata column containing sample identifiers.
#   group_col: Metadata column defining the groups for plotting or summaries.
#   celltype_col: Metadata column containing cell type labels.
#   strip_prefix: Whether to remove sample prefixes from barcodes after splitting or extraction.
# Returns: Prepared sashimi metadata object/table for the next step.
prepare_sashimi_metadata <- function(seurat_obj,
                                     sample_col,
                                     group_col,
                                     celltype_col = NULL,
                                     strip_prefix = TRUE) {
  required_cols <- c(sample_col, group_col)
  if (!is.null(celltype_col) && celltype_col != "") {
    required_cols <- c(required_cols, celltype_col)
  }

  missing_cols <- required_cols[!required_cols %in% colnames(seurat_obj@meta.data)]
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing required metadata columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_barcode") %>%
    mutate(
      sample_id = as.character(.data[[sample_col]]),
      group_id = as.character(.data[[group_col]]),
      plot_barcode = if (strip_prefix) strip_barcode_prefix(cell_barcode) else cell_barcode
    )

  if (!is.null(celltype_col) && celltype_col != "") {
    meta_df$celltype_id <- as.character(meta_df[[celltype_col]])
  }

  meta_df
}

# Brief: Filter sashimi metadata to selected groups, cell types, or samples.
# Args:
#   meta_df: Prepared sashimi metadata table containing cell barcodes, samples, and groups.
#   selected_groups: Group labels to retain for a sashimi comparison.
#   selected_celltypes: Cell type labels to retain for a restricted sashimi comparison.
#   selected_samples: Sample IDs to retain for a restricted sashimi comparison.
# Returns: Sashimi metadata data frame after applying requested group, celltype, and sample filters.
filter_sashimi_metadata <- function(meta_df,
                                    selected_groups = NULL,
                                    selected_celltypes = NULL,
                                    selected_samples = NULL) {
  out <- meta_df

  if (!is.null(selected_groups) && length(selected_groups) > 0) {
    out <- out %>% filter(group_id %in% selected_groups)
  }

  if (!is.null(selected_celltypes) && length(selected_celltypes) > 0) {
    if (!"celltype_id" %in% colnames(out)) {
      stop("selected_celltypes was provided but no celltype metadata was prepared")
    }
    out <- out %>% filter(celltype_id %in% selected_celltypes)
  }

  if (!is.null(selected_samples) && length(selected_samples) > 0) {
    out <- out %>% filter(sample_id %in% selected_samples)
  }

  out
}

# Brief: Extract unique barcodes for every sample-by-group sashimi subset.
# Args:
#   meta_df: Prepared sashimi metadata table containing cell barcodes, samples, and groups.
# Returns: Distinct table of sample_id, group_id, and plot_barcode values for BAM subsetting.
extract_barcodes_by_sample_and_group <- function(meta_df) {
  out <- meta_df %>%
    filter(!is.na(sample_id), sample_id != "", !is.na(group_id), group_id != "") %>%
    distinct(sample_id, group_id, plot_barcode) %>%
    arrange(sample_id, group_id, plot_barcode)

  if (nrow(out) == 0) {
    stop("No barcodes remained after filtering for sashimi plotting")
  }

  out
}

# Brief: Write barcode text files for each sample/group combination.
# Args:
#   barcode_df: Table of sample/group barcodes prepared for sashimi plotting.
#   output_dir: Directory where output files should be written.
# Returns: Manifest data frame for written barcode files.
write_barcode_files <- function(barcode_df, output_dir) {
  ensure_dir(output_dir)

  split_keys <- split(barcode_df, list(barcode_df$sample_id, barcode_df$group_id), drop = TRUE)
  manifest <- lapply(split_keys, function(df) {
    sample_id <- unique(df$sample_id)
    group_id <- unique(df$group_id)
    file_path <- file.path(output_dir, paste0(sample_id, "__", group_id, "_barcodes.txt"))
    writeLines(unique(df$plot_barcode), con = file_path)
    data.frame(
      sample_id = sample_id,
      group_id = group_id,
      barcode_file = file_path,
      n_barcodes = length(unique(df$plot_barcode)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, manifest)
}

# Brief: Build a samtools command for subsetting BAM reads by barcode tag.
# Args:
#   input_bam: Original BAM file to subset.
#   barcode_file: Text file containing barcodes to keep for BAM subsetting.
#   output_bam: Subset BAM file to create or index.
#   barcode_tag: BAM tag used to match cell barcodes.
#   threads: Number of worker threads to use.
# Returns: Shell command string ready to inspect and execute.
prepare_bam_subset_command <- function(input_bam,
                                       barcode_file,
                                       output_bam,
                                       barcode_tag = "CB",
                                       threads = 4) {
  paste(
    "samtools view",
    "-@",
    threads,
    "-b",
    "-D",
    paste0(barcode_tag, ":", shQuote(barcode_file)),
    "-o",
    shQuote(output_bam),
    shQuote(input_bam)
  )
}

# Brief: Build a samtools command for indexing a subset BAM.
# Args:
#   output_bam: Subset BAM file to create or index.
#   threads: Number of worker threads to use.
# Returns: Shell command string ready to inspect and execute.
prepare_bam_index_command <- function(output_bam,
                                      threads = 4) {
  paste(
    "samtools index",
    "-@",
    threads,
    shQuote(output_bam)
  )
}

# Brief: Prepare a manifest describing subset BAM outputs for each sample/group.
# Args:
#   barcode_manifest: Manifest describing barcode files for each sample/group subset.
#   sample_bam_map: Table mapping sample IDs to original BAM paths.
#   subset_bam_dir: Directory where subset BAM files should be written.
# Returns: Prepared subset bam manifest object/table for the next step.
prepare_subset_bam_manifest <- function(barcode_manifest,
                                        sample_bam_map,
                                        subset_bam_dir) {
  ensure_dir(subset_bam_dir)

  required_cols <- c("sample_id", "group_id", "barcode_file")
  missing_cols <- required_cols[!required_cols %in% colnames(barcode_manifest)]
  if (length(missing_cols) > 0) {
    stop(paste("barcode_manifest is missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  if (!all(c("sample_id", "bam_path") %in% colnames(sample_bam_map))) {
    stop("sample_bam_map must contain sample_id and bam_path columns")
  }

  barcode_manifest %>%
    left_join(sample_bam_map, by = "sample_id") %>%
    mutate(
      subset_bam = file.path(subset_bam_dir, paste0(sample_id, "__", group_id, ".bam")),
      sample_name = paste0(sample_id, "__", group_id),
      bam_path = subset_bam,
      group_name = group_id
    )
}

# Brief: Write the three-column ggsashimi BAM input table.
# Args:
#   input_bams_df: Three-column BAM table for ggsashimi input.
#   output_file: Destination path for the generated table or figure.
# Returns: Writes input_bams.tsv and invisibly returns its path.
write_input_bams_tsv <- function(input_bams_df, output_file) {
  required_cols <- c("sample_name", "bam_path", "group_name")
  missing_cols <- required_cols[!required_cols %in% colnames(input_bams_df)]
  if (length(missing_cols) > 0) {
    stop(paste("input_bams_df is missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  ensure_dir(dirname(output_file))
  fwrite(
    input_bams_df[, required_cols, drop = FALSE],
    file = output_file,
    sep = "\t",
    col.names = FALSE
  )

  invisible(output_file)
}

# Brief: Generate default colors for sashimi comparison groups.
# Args:
#   groups: Group labels used to assign colors or plot subsets.
# Returns: Named character vector assigning default colors to group labels.
default_group_palette <- function(groups) {
  base_palette <- c(
    "#1F77B4", "#D62728", "#2CA02C", "#9467BD", "#FF7F0E",
    "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F"
  )

  if (length(groups) > length(base_palette)) {
    base_palette <- rep(base_palette, length.out = length(groups))
  }

  stats::setNames(base_palette[seq_along(groups)], groups)
}

# Brief: Write a ggsashimi-compatible palette file.
# Args:
#   groups: Group labels used to assign colors or plot subsets.
#   output_file: Destination path for the generated table or figure.
#   palette_map: Named vector mapping group labels to colors.
# Returns: Writes palette.txt and invisibly returns its path.
write_palette_txt <- function(groups,
                              output_file,
                              palette_map = NULL) {
  ensure_dir(dirname(output_file))

  if (is.null(palette_map)) {
    palette_map <- default_group_palette(groups)
  }

  missing_groups <- groups[!groups %in% names(palette_map)]
  if (length(missing_groups) > 0) {
    stop(paste("palette_map is missing colors for groups:", paste(missing_groups, collapse = ", ")))
  }

  writeLines(unname(palette_map[groups]), con = output_file)
  invisible(output_file)
}

# Brief: Load transcriptome annotation for gene-region lookup.
# Args:
#   gtf_file: GTF/GFF annotation file used for transcript or gene lookup.
# Returns: Loaded transcriptome gtf object/table from disk or package resources.
read_transcriptome_gtf <- function(gtf_file) {
  rtracklayer::import(gtf_file)
}

# Brief: Resolve the genomic plotting interval for a gene from a GTF file.
# Args:
#   gtf_gr: Imported GTF as a GenomicRanges object.
#   gene_query: Gene symbol or gene ID requested by the user.
#   flank: Number of bases to extend around a resolved gene region.
# Returns: Resolved gene region value or matching records.
resolve_gene_region <- function(gtf_gr,
                                gene_query,
                                flank = 0L) {
  mcols_df <- as.data.frame(mcols(gtf_gr))
  candidate_cols <- intersect(
    c("gene_name", "gene_id", "gene", "gene_symbol"),
    colnames(mcols_df)
  )

  if (length(candidate_cols) == 0) {
    stop("GTF does not contain recognized gene annotation columns such as gene_name or gene_id")
  }

  hits <- rep(FALSE, length(gtf_gr))
  for (col_name in candidate_cols) {
    values <- as.character(mcols_df[[col_name]])
    hits <- hits | (!is.na(values) & values == gene_query)
  }

  gene_gr <- gtf_gr[hits]
  if (length(gene_gr) == 0) {
    stop(paste("No GTF records matched gene query:", gene_query))
  }

  seqname <- as.character(GenomeInfoDb::seqnames(gene_gr)[1])
  start_pos <- max(1L, min(start(gene_gr)) - as.integer(flank))
  end_pos <- max(end(gene_gr)) + as.integer(flank)

  data.frame(
    gene_query = gene_query,
    seqname = seqname,
    start = start_pos,
    end = end_pos,
    region = paste0(seqname, ":", start_pos, "-", end_pos),
    stringsAsFactors = FALSE
  )
}

# Brief: Build the final ggsashimi plotting command.
# Args:
#   input_bams_tsv: Three-column BAM manifest consumed by ggsashimi.
#   region: Genomic interval string used for sashimi plotting.
#   annotation_gtf: Annotation GTF passed to ggsashimi.
#   palette_txt: Palette text file consumed by ggsashimi.
#   min_junction: Minimum junction support threshold for ggsashimi plotting.
#   color_column: Column index in ggsashimi input used for color grouping.
#   overlay_column: Column index in ggsashimi input used for overlay grouping.
#   shrink: Whether ggsashimi should shrink intronic regions in the final plot.
#   alpha: Transparency level passed to ggsashimi for coverage tracks.
#   base_size: Base font size for ggsashimi output.
#   ann_height: Annotation panel height for ggsashimi output.
#   height: Figure height in inches.
#   width: Figure width in inches.
# Returns: Shell command string ready to inspect and execute.
prepare_ggsashimi_command <- function(input_bams_tsv,
                                      region,
                                      annotation_gtf,
                                      palette_txt,
                                      min_junction = 10,
                                      color_column = 3,
                                      overlay_column = 3,
                                      shrink = TRUE,
                                      alpha = 0.25,
                                      base_size = 20,
                                      ann_height = 4,
                                      height = 3,
                                      width = 18) {
  cmd_parts <- c(
    "python /work/scIsoAgent/.pi/scripts/ggsashimi.py",
    "-b", shQuote(input_bams_tsv),
    "-c", shQuote(region),
    "-g", shQuote(annotation_gtf),
    "-M", as.character(min_junction),
    "-C", as.character(color_column),
    "-O", as.character(overlay_column),
    "--alpha", as.character(alpha),
    paste0("--base-size=", base_size),
    paste0("--ann-height=", ann_height),
    paste0("--height=", height),
    paste0("--width=", width),
    "-P", shQuote(palette_txt)
  )

  if (isTRUE(shrink)) {
    cmd_parts <- c(cmd_parts, "--shrink")
  }

  paste(cmd_parts, collapse = " ")
}
