suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DRIMSeq)
  library(DEXSeq)
  library(stageR)
  library(IsoformSwitchAnalyzeR)
  library(rtracklayer)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_dtu_seurat_object <- function(input_rds,
                                   isoform_assay = "isoform") {
  seu.obj <- readRDS(input_rds)

  if (!isoform_assay %in% names(seu.obj@assays)) {
    stop(paste("Isoform assay not found in Seurat object:", isoform_assay))
  }

  DefaultAssay(seu.obj) <- isoform_assay
  seu.obj
}

prepare_flair_dtu_fasta_command <- function(input_fasta,
                                            output_fasta) {
  paste0(
    "awk '\n",
    "/^>/ {\n",
    "    n = split(substr($0, 2), a, \"_\")\n",
    "    if (n > 1) {\n",
    "        header = a[1]\n",
    "        for (i = 2; i < n; i++) {\n",
    "            header = header \"_\" a[i]\n",
    "        }\n",
    "        print \">\" header\n",
    "    } else {\n",
    "        print\n",
    "    }\n",
    "    next\n",
    "}\n",
    "{ print }\n",
    "' ", input_fasta, " > ", output_fasta
  )
}

aggregate_isoform_pseudobulk <- function(seu.obj,
                                         assay = "isoform",
                                         group.by = c("sample_id", "celltype")) {
  pseudo <- AggregateExpression(
    seu.obj,
    assays = assay,
    return.seurat = FALSE,
    group.by = group.by
  )

  as.data.frame(pseudo[[assay]])
}

build_dtu_design_matrix <- function(pseudobulk_df,
                                    split_pattern = "_",
                                    condition_field = 2) {
  sample_data <- data.frame(
    sampleID = colnames(pseudobulk_df),
    stringsAsFactors = FALSE
  )

  split_sample_data <- strsplit(sample_data$sampleID, split_pattern, fixed = TRUE)

  sample_data$condition <- vapply(
    split_sample_data,
    function(x) {
      if (length(x) < condition_field) {
        NA_character_
      } else {
        x[[condition_field]]
      }
    },
    character(1)
  )

  sample_data$condition <- as.factor(sample_data$condition)
  sample_data
}

normalize_transcript_ids_for_dtu <- function(feature_ids) {
  feature_split <- strsplit(feature_ids, "-")
  vapply(
    feature_split,
    function(x) {
      if (length(x) <= 1) {
        x[[1]]
      } else {
        paste(x[-length(x)], collapse = "-")
      }
    },
    character(1)
  )
}

subset_gtf_for_dtu <- function(pseudobulk_df,
                               input_gtf,
                               output_gtf = NULL,
                               normalize_feature_ids = TRUE) {
  tx2gene <- rtracklayer::import(input_gtf)
  cts <- pseudobulk_df

  if (normalize_feature_ids) {
    rownames(cts) <- normalize_transcript_ids_for_dtu(rownames(cts))
  }

  keep_tx <- intersect(tx2gene$transcript_id, rownames(cts))
  subset_tx2gene <- tx2gene[tx2gene$transcript_id %in% keep_tx, ]
  subset_cts <- cts[rownames(cts) %in% keep_tx, , drop = FALSE]

  if (!is.null(output_gtf) && output_gtf != "") {
    ensure_dir(dirname(output_gtf))
    rtracklayer::export(subset_tx2gene, output_gtf, format = "gtf")
  }

  list(
    tx2gene = tx2gene,
    subset_gtf = subset_tx2gene,
    subset_cts = subset_cts
  )
}

initialize_switch_list <- function(subset_cts,
                                   design_matrix,
                                   dtu_gtf,
                                   transcript_fasta,
                                   showProgress = TRUE,
                                   fixStringTieAnnotationProblem = TRUE) {
  importRdata(
    isoformCountMatrix = subset_cts,
    isoformRepExpression = subset_cts,
    designMatrix = design_matrix,
    isoformExonAnnoation = dtu_gtf,
    isoformNtFasta = transcript_fasta,
    showProgress = showProgress,
    fixStringTieAnnotationProblem = fixStringTieAnnotationProblem
  )
}

prefilter_switch_list <- function(aSwitchList,
                                  geneExpressionCutoff = 30,
                                  isoformExpressionCutoff = 20,
                                  IFcutoff = 0.05,
                                  removeSingleIsoformGenes = TRUE,
                                  acceptedGeneBiotype = NULL,
                                  quiet = FALSE) {
  preFilter(
    aSwitchList,
    geneExpressionCutoff = geneExpressionCutoff,
    isoformExpressionCutoff = isoformExpressionCutoff,
    IFcutoff = IFcutoff,
    removeSingleIsoformGenes = removeSingleIsoformGenes,
    acceptedGeneBiotype = acceptedGeneBiotype,
    quiet = quiet
  )
}

add_orf_and_gene_symbols <- function(aSwitchList_filt,
                                     resource_table_path,
                                     orfMethod = "longest",
                                     showProgress = TRUE) {
  aSwitchList_ORF_filt <- analyzeORF(
    switchAnalyzeRlist = aSwitchList_filt,
    orfMethod = orfMethod,
    showProgress = showProgress
  )

  resource_table <- read.csv(resource_table_path, header = TRUE, stringsAsFactors = FALSE)

  if (!all(c("gene_id", "gene_symbol") %in% colnames(resource_table))) {
    stop("resource table must contain gene_id and gene_symbol columns")
  }

  isoform_features <- aSwitchList_ORF_filt[["isoformFeatures"]]

  # Prefer transcript-level mapping to avoid ID-system mismatch
  # (e.g., switch list gene_id may be symbols while resource gene_id is ENSG).
  mapped_gene_symbol <- rep(NA_character_, nrow(isoform_features))

  if ("transcript_id" %in% colnames(resource_table) && "isoform_id" %in% colnames(isoform_features)) {
    rt <- resource_table %>%
      dplyr::transmute(
        transcript_id = as.character(transcript_id),
        gene_id = as.character(gene_id),
        gene_symbol = as.character(gene_symbol)
      ) %>%
      dplyr::filter(!is.na(transcript_id), transcript_id != "") %>%
      dplyr::distinct(transcript_id, .keep_all = TRUE)

    idx <- match(as.character(isoform_features[["isoform_id"]]), rt$transcript_id)
    mapped_gene_symbol <- rt$gene_symbol[idx]

    # Fallback: strip transcript version suffix if needed
    miss <- is.na(mapped_gene_symbol) | mapped_gene_symbol == ""
    if (any(miss)) {
      idx2 <- match(
        sub("\\.[0-9]+$", "", as.character(isoform_features[["isoform_id"]])[miss]),
        sub("\\.[0-9]+$", "", rt$transcript_id)
      )
      mapped_gene_symbol[miss] <- rt$gene_symbol[idx2]
    }
  }

  # Secondary fallback by gene_id matching when still missing
  miss2 <- is.na(mapped_gene_symbol) | mapped_gene_symbol == ""
  if (any(miss2) && "gene_id" %in% colnames(isoform_features)) {
    mapped_gene_symbol[miss2] <- resource_table$gene_symbol[
      match(as.character(isoform_features[["gene_id"]])[miss2], as.character(resource_table$gene_id))
    ]
  }

  isoform_features[["gene_name"]] <- mapped_gene_symbol
  aSwitchList_ORF_filt[["isoformFeatures"]] <- isoform_features

  aSwitchList_ORF_filt
}

run_dtu_test <- function(aSwitchList_ORF_filt,
                         reduceToSwitchingGenes = FALSE,
                         reduceFurtherToGenesWithConsequencePotential = FALSE,
                         alpha = 0.01,
                         dIFcutoff = 0.25,
                         onlySigIsoforms = FALSE) {
  isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchList_ORF_filt,
    reduceToSwitchingGenes = reduceToSwitchingGenes,
    reduceFurtherToGenesWithConsequencePotential = reduceFurtherToGenesWithConsequencePotential,
    alpha = alpha,
    dIFcutoff = dIFcutoff,
    onlySigIsoforms = onlySigIsoforms
  )
}

export_dtu_results <- function(SwitchListAnalyzed,
                               output_csv,
                               top_n = 20) {
  ensure_dir(dirname(output_csv))

  iso_feat <- as.data.frame(SwitchListAnalyzed$isoformFeatures)
  write.csv(iso_feat, output_csv, row.names = FALSE)

  list(
    switch_summary = extractSwitchSummary(SwitchListAnalyzed),
    top_switches = extractTopSwitches(SwitchListAnalyzed, n = top_n),
    isoform_features = iso_feat
  )
}

plot_dtu_switch <- function(SwitchListAnalyzed,
                            gene,
                            condition1,
                            condition2,
                            output_file = NULL,
                            width = 14,
                            height = 8,
                            localTheme = theme_bw(base_size = 12)) {
  if (!is.null(output_file) && output_file != "") {
    ensure_dir(dirname(output_file))
    pdf(output_file, width = width, height = height, bg = "white")
    print(
      switchPlot(
        SwitchListAnalyzed,
        gene = gene,
        condition1 = condition1,
        condition2 = condition2,
        localTheme = localTheme
      )
    )
    dev.off()
    return(invisible(NULL))
  }

  switchPlot(
    SwitchListAnalyzed,
    gene = gene,
    condition1 = condition1,
    condition2 = condition2,
    localTheme = localTheme
  )
}
