suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

load_annotation_libraries <- function(extra_packages = c("HGNChelper", "ggraph", "igraph", "tidyverse", "data.tree")) {
  invisible(lapply(extra_packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
    }
  }))
}

load_sctype_functions <- function() {
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_plot_pdf <- function(plot_obj, file_path, width = 8, height = 8) {
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

read_annotation_object <- function(input_rds) {
  readRDS(input_rds)
}

plot_cluster_umap <- function(seurat_obj,
                              group_by,
                              output_file = NULL,
                              reduction = "umap",
                              title = "Cluster UMAP") {
  p <- CellDimPlot(
    srt = seurat_obj,
    group.by = group_by,
    reduction = reduction,
    theme_use = "theme_blank",
    title = title
  )

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    save_plot_pdf(p, output_file, width = 8, height = 8)
  }

  p
}

compute_cluster_markers <- function(seurat_obj,
                                    assay = "RNA",
                                    logfc.threshold = 0.5,
                                    min.pct = 0.20,
                                    only.pos = TRUE,
                                    p_adj_cutoff = 0.05,
                                    output_file = NULL) {
  markers <- FindAllMarkers(
    seurat_obj,
    assay = assay,
    do.print = TRUE,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    only.pos = only.pos
  ) %>%
    dplyr::filter(p_val_adj < p_adj_cutoff)

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    write.csv(markers, output_file, row.names = FALSE)
  }

  markers
}

perform_sctype_annotation <- function(seurat_obj,
                                      db_,
                                      tissue,
                                      gs_removal_list = character(0),
                                      metadata_col_prefix = "sctype_db",
                                      figure_prefix = "annotation",
                                      cluster_res = "RNA_snn_res.0.9",
                                      reduction = "umap",
                                      output_pdf = NULL) {
  if (!exists("gene_sets_prepare", mode = "function") || !exists("sctype_score", mode = "function")) {
    load_sctype_functions()
  }

  gs_list <- gene_sets_prepare(db_, tissue)

  for (gs in gs_removal_list) {
    if (gs %in% names(gs_list[["gs_positive"]])) {
      gs_list[["gs_positive"]][[gs]] <- NULL
    }
  }

  es.max <- sctype_score(
    scRNAseqData = seurat_obj[["RNA"]]@scale.data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )

  Idents(seurat_obj) <- cluster_res

  cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data[[cluster_res]]), function(cl) {
    es.max.cl <- sort(
      rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_res]] == cl, , drop = FALSE])]),
      decreasing = TRUE
    )
    head(
      data.frame(
        cluster = cl,
        type = names(es.max.cl),
        scores = es.max.cl,
        ncells = sum(seurat_obj@meta.data[[cluster_res]] == cl)
      ),
      10
    )
  }))

  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores) %>%
    as.data.frame()

  sctype_scores$scores <- as.numeric(sctype_scores$scores)
  sctype_scores$type[sctype_scores$scores < sctype_scores$ncells / 4] <- "Unknown"

  seurat_obj@meta.data[[metadata_col_prefix]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, , drop = FALSE]
    seurat_obj@meta.data[[metadata_col_prefix]][seurat_obj@meta.data[[cluster_res]] == j] <- as.character(cl_type$type[1])
  }

  p <- CellDimPlot(
    srt = seurat_obj,
    group.by = metadata_col_prefix
    ,
    reduction = reduction,
    theme_use = "theme_blank",
    title = figure_prefix
  )

  if (!is.null(output_pdf) && output_pdf != "") {
    ensure_dir(dirname(output_pdf))
    save_plot_pdf(p, output_pdf, width = 8, height = 8)
  }

  list(
    seurat_obj = seurat_obj,
    sctype_scores = sctype_scores,
    plot = p
  )
}

rename_clusters_manually <- function(seurat_obj,
                                     cluster_to_celltype,
                                     cluster_col = NULL,
                                     celltype_col = "celltype") {
  if (!is.null(cluster_col) && cluster_col != "") {
    Idents(seurat_obj) <- cluster_col
  }

  seurat_obj <- RenameIdents(seurat_obj, !!!cluster_to_celltype)
  seurat_obj[[celltype_col]] <- Idents(seurat_obj)
  seurat_obj
}

plot_celltype_umap <- function(seurat_obj,
                               group_by = "celltype",
                               reduction = "umap",
                               output_file = NULL,
                               title = "Celltype UMAP") {
  p <- CellDimPlot(
    srt = seurat_obj,
    group.by = group_by,
    reduction = reduction,
    theme_use = "theme_blank",
    title = title
  )

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    save_plot_pdf(p, output_file, width = 8, height = 8)
  }

  p
}

run_celltype_de <- function(seurat_obj,
                            group_by = "celltype",
                            assay = "RNA",
                            fc.threshold = 1,
                            only.pos = FALSE,
                            output_file = NULL) {
  seurat_obj[[group_by]] <- seurat_obj@meta.data[[group_by]]
  seurat_obj <- RunDEtest(
    srt = seurat_obj,
    group_by = group_by,
    fc.threshold = fc.threshold,
    only.pos = only.pos,
    assay = assay
  )

  de_table <- seurat_obj@tools$DEtest_celltype$AllMarkers_wilcox

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    write.csv(de_table, output_file, row.names = FALSE)
  }

  list(
    seurat_obj = seurat_obj,
    de_table = de_table
  )
}

plot_volcano <- function(seurat_obj,
                         group_by = "celltype",
                         output_file = NULL,
                         width = 14,
                         height = 16) {
  p_list <- VolcanoPlot(srt = seurat_obj, group_by = group_by, combine = FALSE, pt.size = 0.3)
  p_list2 <- lapply(seq_along(p_list), function(i) {
    p_list[[i]] +
      theme(legend.position = "none") +
      labs(title = names(p_list)[i])
  })
  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    save_plot_pdf(cowplot::plot_grid(plotlist = p_list2, ncol = 2), output_file, width = width, height = height)
  }
  cowplot::plot_grid(plotlist = p_list2, ncol = 2)
}

plot_celltype_heatmap <- function(seurat_obj,
                                  de_table,
                                  group_by = "celltype",
                                  species = "Homo_sapiens",
                                  assay = "RNA",
                                  db = c("TF", "CSPA"),
                                  fc_cutoff = 1,
                                  p_adj_cutoff = 0.05,
                                  anno_terms = assay == "RNA",
                                  output_file = NULL,
                                  height = 5,
                                  width = 4) {
  if (!(assay %in% names(seurat_obj@assays))) {
    stop(sprintf("Assay '%s' not found in Seurat object", assay))
  }
  DefaultAssay(seurat_obj) <- assay

  de_filtered <- de_table[with(de_table, avg_log2FC > fc_cutoff & p_val_adj < p_adj_cutoff), , drop = FALSE]
  if (nrow(de_filtered) == 0) {
    stop("No differential features passed the requested thresholds for heatmap plotting")
  }

  heatmap_features <- if ("gene" %in% colnames(de_filtered)) de_filtered$gene else rownames(de_filtered)
  heatmap_features <- unique(as.character(heatmap_features))
  heatmap_features <- heatmap_features[heatmap_features %in% rownames(seurat_obj[[assay]])]
  if (length(heatmap_features) == 0) {
    stop(sprintf("No requested features were found in assay '%s'", assay))
  }

  feature_split <- NULL
  if ("group1" %in% colnames(de_filtered)) {
    idx <- match(heatmap_features, as.character(if ("gene" %in% colnames(de_filtered)) de_filtered$gene else rownames(de_filtered)))
    feature_split <- as.character(de_filtered$group1[idx])
  }

  db_use <- intersect(db, c("TF", "CSPA"))
  if (length(db_use) > 0) {
    seurat_obj <- AnnotateFeatures(seurat_obj, species = species, db = db_use)
  }

  feature_meta_cols <- colnames(seurat_obj[[assay]]@meta.features)
  feature_annotation_use <- NULL
  feature_annotation_pal_use <- NULL
  if (assay == "RNA") {
    feature_annotation_use <- intersect(c("TF", "CSPA"), feature_meta_cols)
    if (length(feature_annotation_use) > 0) {
      pal_map <- list(TF = c("gold", "steelblue"), CSPA = c("forestgreen"))
      feature_annotation_pal_use <- unname(pal_map[feature_annotation_use])
    }
  }

  ht <- FeatureHeatmap(
    srt = seurat_obj,
    group.by = group_by,
    features = heatmap_features,
    feature_split = feature_split,
    species = species,
    db = if (anno_terms) c("GO_BP", "KEGG", "WikiPathway") else NULL,
    anno_terms = anno_terms,
    feature_annotation = feature_annotation_use,
    feature_annotation_palcolor = feature_annotation_pal_use,
    height = height,
    width = width
  )

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    pdf(output_file, width = width + 6, height = height + 4, bg = "white")
    if (is.list(ht) && !is.null(ht$plot)) {
      print(ht$plot)
    } else {
      print(ht)
    }
    dev.off()
  }

  list(
    seurat_obj = seurat_obj,
    heatmap = ht
  )
}

run_celltype_enrichment <- function(seurat_obj,
                                    group_by = "celltype",
                                    db = "GO_BP",
                                    species = "Homo_sapiens",
                                    de_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05") {
  RunEnrichment(
    srt = seurat_obj,
    group_by = group_by,
    db = db,
    species = species,
    DE_threshold = de_threshold
  )
}

plot_enrichment <- function(seurat_obj,
                            group_by = "celltype",
                            group_use,
                            plot_type = "bar",
                            word_type = NULL,
                            output_file = NULL,
                            width = 10,
                            height = 8) {
  args <- list(
    srt = seurat_obj,
    group_by = group_by,
    group_use = group_use,
    plot_type = plot_type
  )
  if (!is.null(word_type)) {
    args$word_type <- word_type
  }

  p <- do.call(EnrichmentPlot, args)

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    save_plot_pdf(p, output_file, width = width, height = height)
  }

  p
}

plot_cell_statistics <- function(seurat_obj,
                                 stat.by,
                                 group.by,
                                 label = TRUE,
                                 output_file = NULL,
                                 width = 10,
                                 height = 8) {
  p <- CellStatPlot(
    srt = seurat_obj,
    stat.by = stat.by,
    group.by = group.by,
    label = label
  )

  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    save_plot_pdf(p, output_file, width = width, height = height)
  }

  p
}
