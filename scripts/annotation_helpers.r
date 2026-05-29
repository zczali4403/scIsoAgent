# ------------------------------------------------------------------------------
# Script: annotation_helpers.r
# Purpose: Helper utilities for downstream cell-type annotation and marker-based summaries.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# Brief: Load optional annotation-related R packages when they are installed.
# Args:
#   extra_packages: Optional package names to load when available.
# Returns: Loaded annotation libraries object/table from disk or package resources.
load_annotation_libraries <- function(extra_packages = c("HGNChelper", "ggraph", "igraph", "tidyverse", "data.tree")) {
  invisible(lapply(extra_packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
    }
  }))
}

# Brief: Load ScType helper functions from their upstream source URLs.
# Args: None.
# Returns: Loaded sctype functions object/table from disk or package resources.
load_sctype_functions <- function() {
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
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
save_plot_pdf <- function(plot_obj, file_path, width = 8, height = 8) {
  pdf(file = file_path, width = width, height = height, bg = "white")
  print(plot_obj)
  dev.off()
}

# Brief: Load the Seurat object to annotate.
# Args:
#   input_rds: Path to the Seurat RDS object used as input.
# Returns: Loaded Seurat object for annotation.
read_annotation_object <- function(input_rds) {
  readRDS(input_rds)
}

# Brief: Plot clusters on an embedding before or during annotation.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   output_file: Destination path for the generated table or figure.
#   reduction: Dimensional reduction name used for embedding plots.
#   title: Plot title text.
# Returns: Plot object for cluster umap; writes a PDF when an output path is provided.
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

# Brief: Find cluster marker features for a selected assay and threshold set.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   logfc.threshold: Minimum log fold-change threshold for marker detection.
#   min.pct: Minimum fraction of cells expressing a marker feature.
#   only.pos: Whether marker detection should return only positive markers.
#   p_adj_cutoff: Adjusted p-value cutoff for filtering marker or DE results.
#   output_file: Destination path for the generated table or figure.
# Returns: Computed cluster markers table or vector.
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

# Brief: Run ScType scoring and add predicted cell-type labels to Seurat metadata.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   db_: ScType marker database identifier or path.
#   tissue: Tissue name used to select ScType marker sets.
#   gs_removal_list: ScType gene-set names to remove before scoring.
#   metadata_col_prefix: Metadata column name used to store ScType labels.
#   figure_prefix: Text prefix used in annotation plot titles or filenames.
#   cluster_res: Clustering resolution used by Seurat FindClusters.
#   reduction: Dimensional reduction name used for embedding plots.
#   output_pdf: PDF path for the aggregated Orthrus report.
# Returns: List containing updated Seurat object, ScType score table, and annotation UMAP plot.
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

# Brief: Apply user-provided cluster-to-celltype labels to Seurat identities and metadata.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   cluster_to_celltype: Named mapping from cluster IDs to final cell type labels.
#   cluster_col: Metadata column containing cluster identities.
#   celltype_col: Metadata column containing cell type labels.
# Returns: Seurat object with identities renamed and the celltype metadata column updated.
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

# Brief: Plot annotated cell types on a Seurat embedding.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   reduction: Dimensional reduction name used for embedding plots.
#   output_file: Destination path for the generated table or figure.
#   title: Plot title text.
# Returns: Plot object for celltype umap; writes a PDF when an output path is provided.
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

# Brief: Run celltype-level differential expression or differential isoform testing with SCP.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   fc.threshold: Fold-change threshold for SCP differential testing.
#   only.pos: Whether marker detection should return only positive markers.
#   output_file: Destination path for the generated table or figure.
# Returns: List containing the updated Seurat object and the celltype DE marker table.
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

# Brief: Create per-celltype volcano plots from SCP differential-test results.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   output_file: Destination path for the generated table or figure.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Plot object for volcano; writes a PDF when an output path is provided.
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

# Brief: Generate feature heatmaps for celltype-level marker or DE features.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   de_table: Differential expression/feature table used to select plot features.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   species: Species name used by SCP annotation/enrichment helpers.
#   assay: Seurat assay name or SingleCellExperiment assay name to use.
#   db: Database name(s) used for annotation, enrichment, or feature annotation.
#   fc_cutoff: Fold-change cutoff used to select features for plotting.
#   p_adj_cutoff: Adjusted p-value cutoff for filtering marker or DE results.
#   anno_terms: Whether to add pathway or ontology annotations to a heatmap.
#   output_file: Destination path for the generated table or figure.
#   height: Figure height in inches.
#   width: Figure width in inches.
# Returns: List containing the updated Seurat object and heatmap object.
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

# Brief: Run pathway or gene-set enrichment analysis for celltype DE results.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   db: Database name(s) used for annotation, enrichment, or feature annotation.
#   species: Species name used by SCP annotation/enrichment helpers.
#   de_threshold: Expression used by SCP to filter DE features for enrichment.
# Returns: Result object/table from the celltype enrichment step.
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

# Brief: Plot enrichment results for selected celltype groups.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   group_by: Metadata column(s) used to aggregate, compare, or plot cell groups.
#   group_use: Selected group(s) to display in an enrichment plot.
#   plot_type: Enrichment plot style to produce.
#   word_type: Wordcloud term type used by SCP enrichment plotting.
#   output_file: Destination path for the generated table or figure.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Plot object for enrichment; writes a PDF when an output path is provided.
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

# Brief: Plot cell composition summaries across metadata groupings.
# Args:
#   seurat_obj: Seurat object being read, modified, plotted, or summarized.
#   stat.by: Metadata variable used for cell composition statistics.
#   group.by: Metadata variable used for grouping cells in Seurat/SCP calls.
#   label: Whether labels should be shown on a composition plot.
#   output_file: Destination path for the generated table or figure.
#   width: Figure width in inches.
#   height: Figure height in inches.
# Returns: Plot object for cell statistics; writes a PDF when an output path is provided.
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
