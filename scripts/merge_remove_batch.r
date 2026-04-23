suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(data.table)
  library(ggplot2)
})

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

require_args <- function(parsed, required_keys) {
  missing_keys <- required_keys[!required_keys %in% names(parsed)]
  if (length(missing_keys) > 0) {
    stop(
      paste(
        "Missing required arguments:",
        paste(missing_keys, collapse = ", "),
        "\nRequired usage:",
        paste(
          "Rscript scripts/merge_remove_batch.r",
          "--stage=elbow|preview|harmony-elbow|finalize",
          "--input_rds_files=/path/s1_object.rds,/path/s2_object.rds",
          "--output_dir=/path/to/output_dir",
          "--output_prefix=merged_samples",
          "[--sample_names=sample1,sample2]",
          "[--object_label=gene|isoform|custom_label]",
          "[--remove_batch=false]",
          "[--batch_method=harmony]",
          "[--npc=30]",
          "[--harmony_npc=20]",
          "[--cluster_res=0.5]"
        )
      )
    )
  }
}

split_csv_arg <- function(value) {
  if (is.null(value) || value == "") {
    return(character(0))
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

parse_bool_arg <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(default)
  }
  tolower(value) %in% c("true", "1", "yes", "y")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

infer_sample_names <- function(paths) {
  sample_names <- basename(paths)
  sample_names <- sub("\\.rds$", "", sample_names, ignore.case = TRUE)
  sample_names <- sub("_isoform_subset$", "", sample_names, ignore.case = TRUE)
  sample_names <- sub("_gene_singlets$", "", sample_names, ignore.case = TRUE)
  sample_names <- sub("_singlets$", "", sample_names, ignore.case = TRUE)
  sample_names <- sub("^removed_empty_", "", sample_names, ignore.case = TRUE)
  sample_names <- gsub("[^A-Za-z0-9._-]", "_", sample_names)

  if (any(sample_names == "")) {
    stop("Failed to infer sample names from input_rds_files")
  }
  if (anyDuplicated(sample_names)) {
    stop("Inferred sample names are not unique; please provide --sample_names explicitly")
  }

  sample_names
}

validate_input_lengths <- function(input_rds_files, sample_names) {
  if (length(input_rds_files) == 0) {
    stop("At least one input RDS file must be provided")
  }
  if (length(input_rds_files) != length(sample_names)) {
    stop("input_rds_files and sample_names must have the same number of entries")
  }
}

load_seurat_objects <- function(input_rds_files, sample_names) {
  objects <- vector("list", length(input_rds_files))

  for (i in seq_along(input_rds_files)) {
    obj <- readRDS(input_rds_files[[i]])
    obj$sample_id <- sample_names[[i]]
    obj$orig.ident <- sample_names[[i]]
    objects[[i]] <- obj
  }

  objects
}

merge_objects <- function(objects, sample_names) {
  if (length(objects) == 1) {
    return(objects[[1]])
  }

  merge(
    x = objects[[1]],
    y = objects[-1],
    add.cell.ids = sample_names,
    project = "merged_samples"
  )
}

prepare_merged_object <- function(merged_obj, npc) {
  merged_obj <- merged_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npc, verbose = FALSE)
  merged_obj
}

build_elbow_plot <- function(merged_obj) {
  ElbowPlot(merged_obj, ndims = 50) +
    ggtitle("Elbow Plot") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

build_harmony_elbow_plot <- function(merged_obj) {
  if (!"harmony" %in% names(merged_obj@reductions)) {
    stop("Harmony reduction not found; run Harmony before building a harmony elbow plot")
  }

  harmony_embeddings <- Embeddings(merged_obj, reduction = "harmony")
  stdev <- apply(harmony_embeddings, 2, stats::sd)
  elbow_df <- data.frame(
    dims = seq_along(stdev),
    stdev = stdev
  )

  ggplot(elbow_df, aes(x = dims, y = stdev)) +
    geom_point(color = "#2563EB", size = 1.5) +
    geom_line(color = "#2563EB", linewidth = 0.5) +
    labs(
      title = "Harmony Elbow Plot",
      x = "Harmony Dimension",
      y = "Standard Deviation"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

run_batch_correction <- function(merged_obj, batch_method, npc) {
  if (batch_method != "harmony") {
    stop("Unsupported --batch_method. Currently supported: harmony")
  }

  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Harmony batch correction requested, but the harmony package is not installed")
  }

  merged_obj <- harmony::RunHarmony(
    object = merged_obj,
    group.by.vars = "sample_id",
    reduction.use = "pca",
    dims.use = 1:npc
  )

  merged_obj
}

run_graph_and_umap <- function(merged_obj, reduction_name, dims_to_use, cluster_res) {
  merged_obj <- FindNeighbors(merged_obj, reduction = reduction_name, dims = 1:dims_to_use)
  merged_obj <- FindClusters(merged_obj, resolution = cluster_res)
  merged_obj <- RunUMAP(merged_obj, reduction = reduction_name, dims = 1:dims_to_use)
  merged_obj
}

save_dimplot <- function(merged_obj, group_by, file_path, title_text) {
  p <- CellDimPlot(
    srt = merged_obj,
    group.by = group_by,
    reduction = "UMAP",
    theme_use = "theme_blank",
    title = title_text
  )

  pdf(file_path, width = 7, height = 6, bg = "white")
  print(p)
  dev.off()
}

write_run_summary <- function(output_dir,
                              output_prefix,
                              object_label,
                              remove_batch,
                              batch_method,
                              npc,
                              cluster_res,
                              sample_names,
                              merged_obj,
                              harmony_npc = NA) {
  run_summary <- data.frame(
    parameter = c("object_label", "remove_batch", "batch_method", "npc", "harmony_npc", "cluster_res", "n_samples", "merged_cells"),
    value = c(object_label, remove_batch, batch_method, npc, harmony_npc, cluster_res, length(sample_names), ncol(merged_obj)),
    stringsAsFactors = FALSE
  )
  fwrite(run_summary, file.path(output_dir, paste0(output_prefix, "_run_summary.tsv")), sep = "\t")
}

build_merged_object <- function(input_rds_files,
                                sample_names,
                                npc,
                                cluster_res = NULL,
                                remove_batch,
                                batch_method,
                                run_graph = TRUE,
                                graph_dims = NULL) {
  objects <- load_seurat_objects(input_rds_files, sample_names)
  merged_obj <- merge_objects(objects, sample_names)
  merged_obj <- prepare_merged_object(merged_obj, npc)

  reduction_name <- "pca"
  if (remove_batch) {
    merged_obj <- run_batch_correction(merged_obj, batch_method, npc)
    reduction_name <- if (batch_method == "harmony") "harmony" else reduction_name
  }

  if (run_graph) {
    dims_to_use <- if (is.null(graph_dims)) npc else graph_dims
    merged_obj <- run_graph_and_umap(merged_obj, reduction_name, dims_to_use, cluster_res)
  }

  list(
    objects = objects,
    merged_obj = merged_obj,
    reduction_name = reduction_name
  )
}

write_sample_summary <- function(output_dir, output_prefix, sample_names, input_rds_files, objects) {
  summary_df <- data.frame(
    sample_name = sample_names,
    input_rds_file = input_rds_files,
    cells = vapply(objects, ncol, integer(1)),
    features = vapply(objects, nrow, integer(1)),
    stringsAsFactors = FALSE
  )
  fwrite(summary_df, file.path(output_dir, paste0(output_prefix, "_sample_summary.tsv")), sep = "\t")
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_args(args, c("stage", "input_rds_files", "output_dir", "output_prefix"))

  stage <- args$stage
  input_rds_files <- split_csv_arg(args$input_rds_files)
  sample_names <- if ("sample_names" %in% names(args)) {
    split_csv_arg(args$sample_names)
  } else {
    infer_sample_names(input_rds_files)
  }
  validate_input_lengths(input_rds_files, sample_names)

  output_dir <- args$output_dir
  output_prefix <- args$output_prefix
  remove_batch <- parse_bool_arg(args$remove_batch, default = FALSE)
  batch_method <- if ("batch_method" %in% names(args)) args$batch_method else "harmony"
  object_label <- if ("object_label" %in% names(args)) args$object_label else "seurat_object"
  npc <- if ("npc" %in% names(args)) as.integer(args$npc) else NA_integer_
  harmony_npc <- if ("harmony_npc" %in% names(args)) as.integer(args$harmony_npc) else NA_integer_
  cluster_res <- if ("cluster_res" %in% names(args)) as.numeric(args$cluster_res) else 0.5

  ensure_dir(output_dir)

  if (!stage %in% c("elbow", "preview", "harmony-elbow", "finalize")) {
    stop("Unsupported --stage. Use one of: elbow, preview, harmony-elbow, finalize")
  }

  objects <- load_seurat_objects(input_rds_files, sample_names)
  merged_obj <- merge_objects(objects, sample_names)
  merged_obj <- prepare_merged_object(merged_obj, npc = 50L)
  write_sample_summary(output_dir, output_prefix, sample_names, input_rds_files, objects)

  if (stage == "elbow") {
    elbow_plot <- build_elbow_plot(merged_obj)
    pdf(file.path(output_dir, paste0(output_prefix, "_elbowplot.pdf")), width = 7, height = 6, bg = "white")
    print(elbow_plot)
    dev.off()
    write_run_summary(
      output_dir = output_dir,
      output_prefix = output_prefix,
      object_label = object_label,
      remove_batch = FALSE,
      batch_method = "pending",
      npc = "pending",
      harmony_npc = "pending",
      cluster_res = cluster_res,
      sample_names = sample_names,
      merged_obj = merged_obj
    )
    return(invisible(NULL))
  }

  if (is.na(npc)) {
    stop("--npc is required for preview, harmony-elbow, and finalize stages after the user reviews the PCA elbow plot")
  }

  if (stage == "harmony-elbow") {
    built <- build_merged_object(
      input_rds_files = input_rds_files,
      sample_names = sample_names,
      npc = npc,
      remove_batch = TRUE,
      batch_method = batch_method,
      run_graph = FALSE
    )
    merged_obj <- built$merged_obj
    harmony_elbow_plot <- build_harmony_elbow_plot(merged_obj)
    pdf(file.path(output_dir, paste0(output_prefix, "_harmony_elbowplot.pdf")), width = 7, height = 6, bg = "white")
    print(harmony_elbow_plot)
    dev.off()
    write_run_summary(
      output_dir = output_dir,
      output_prefix = output_prefix,
      object_label = object_label,
      remove_batch = TRUE,
      batch_method = paste0(batch_method, "_pending_dims"),
      npc = paste0("pca_npc=", npc, "; harmony_npc=pending"),
      harmony_npc = "pending",
      cluster_res = cluster_res,
      sample_names = sample_names,
      merged_obj = merged_obj
    )
    return(invisible(NULL))
  }

  built <- build_merged_object(
    input_rds_files = input_rds_files,
    sample_names = sample_names,
    npc = npc,
    cluster_res = cluster_res,
    remove_batch = if (stage == "preview") FALSE else remove_batch,
    batch_method = batch_method,
    run_graph = TRUE,
    graph_dims = if (stage == "finalize" && remove_batch) {
      if (is.na(harmony_npc)) stop("--harmony_npc is required for finalize when --remove_batch=true after reviewing the harmony elbow plot")
      harmony_npc
    } else {
      NULL
    }
  )
  merged_obj <- built$merged_obj

  if (stage == "preview") {
    write_run_summary(
      output_dir = output_dir,
      output_prefix = output_prefix,
      object_label = object_label,
      remove_batch = FALSE,
      batch_method = "none",
      npc = npc,
      harmony_npc = "not_used",
      cluster_res = cluster_res,
      sample_names = sample_names,
      merged_obj = merged_obj
    )
    save_dimplot(
      merged_obj,
      "sample_id",
      file.path(output_dir, paste0(output_prefix, "_preview_by_sample_umap.pdf")),
      paste0("Preview ", object_label, " UMAP by Sample (No Batch Correction)")
    )
    save_dimplot(
      merged_obj,
      "seurat_clusters",
      file.path(output_dir, paste0(output_prefix, "_preview_by_cluster_umap.pdf")),
      paste0("Preview ", object_label, " UMAP by Cluster (No Batch Correction)")
    )
    return(invisible(NULL))
  }

  write_run_summary(
    output_dir = output_dir,
    output_prefix = output_prefix,
    object_label = object_label,
    remove_batch = remove_batch,
    batch_method = if (remove_batch) batch_method else "none",
    npc = npc,
    harmony_npc = if (remove_batch) harmony_npc else "not_used",
    cluster_res = cluster_res,
    sample_names = sample_names,
    merged_obj = merged_obj
  )
  saveRDS(merged_obj, file.path(output_dir, paste0(output_prefix, ".rds")))
  save_dimplot(
    merged_obj,
    "sample_id",
    file.path(output_dir, paste0(output_prefix, "_by_sample_umap.pdf")),
    if (remove_batch) {
      paste0("Merged ", object_label, " UMAP by Sample (Batch Corrected)")
    } else {
      paste0("Merged ", object_label, " UMAP by Sample")
    }
  )
  save_dimplot(
    merged_obj,
    "seurat_clusters",
    file.path(output_dir, paste0(output_prefix, "_by_cluster_umap.pdf")),
    if (remove_batch) {
      paste0("Merged ", object_label, " UMAP by Cluster (Batch Corrected)")
    } else {
      paste0("Merged ", object_label, " UMAP by Cluster")
    }
  )
}

main()
