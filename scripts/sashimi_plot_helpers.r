suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(rtracklayer)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_sashimi_seurat_object <- function(input_rds) {
  readRDS(input_rds)
}

list_metadata_columns <- function(seurat_obj) {
  colnames(seurat_obj@meta.data)
}

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

strip_barcode_prefix <- function(barcodes) {
  sub("^[^_]+_", "", barcodes)
}

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

prepare_bam_index_command <- function(output_bam,
                                      threads = 4) {
  paste(
    "samtools index",
    "-@",
    threads,
    shQuote(output_bam)
  )
}

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

read_transcriptome_gtf <- function(gtf_file) {
  rtracklayer::import(gtf_file)
}

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
