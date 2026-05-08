---
name: annotation
description: Performs downstream cell type annotation for the project in one persistent tmux R session, including cluster UMAP review, cluster-level marker discovery for gene and isoform assays, optional automatic annotation with ScType, optional manual cluster renaming, celltype-level differential analysis, heatmaps, enrichment, and cell composition summaries.
---

# Annotation

Use this skill when the user wants downstream cell type annotation after a merged Seurat object has already been prepared.

This skill is for annotation and annotation-driven downstream interpretation only. It does not cover upstream quantification, empty-droplet removal, QC, object merging, or trajectory analysis.

## Use This Skill When

- the user wants to annotate clusters into cell types
- the user wants automatic annotation
- the user wants manual annotation cluster by cluster
- the user wants cluster markers or celltype markers to support annotation
- the user wants to inspect gene or isoform feature plots during annotation

## Mandatory Session Rule

All analysis in this skill must be performed in one persistent interactive `tmux` R session.

Do not split this workflow across repeated `Rscript` calls.
Do not exit the R session after one substep.
Keep the same `tmux` R terminal alive until the user explicitly says to exit annotation.

## Project Source Pattern

This skill is primarily an interactive workflow skill.

It relies on standard R code inside the persistent session, including:

- `Seurat`
- `SCP`
- `HGNChelper`
- `ggraph`
- `igraph`
- `tidyverse`
- `data.tree`

For automatic annotation, it also uses:

- `gene_sets_prepare.R`
- `sctype_score_.R`

Primary local helper script:

- `/work/.pi/scripts/annotation_helpers.r`

loaded from:

- `https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R`
- `https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R`

At the start of the annotation session, prefer:

```r
source("/work/.pi/scripts/annotation_helpers.r")
```

Then call helper functions from that script instead of repeatedly rewriting the same annotation code in the terminal.

## Entry Procedure

At the start of an annotation request:

1. Use the `downstream` conda environment for the whole annotation workflow.
2. Start one persistent `tmux` R session and keep it open for the whole annotation workflow.
3. Resolve the merged Seurat RDS path to annotate.
4. In that session, `source("/work/.pi/scripts/annotation_helpers.r")`.
5. Read the Seurat object once in that session.
6. Derive the output directory for annotation artifacts from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
7. Use `13_annotation` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
8. Confirm which clustering metadata column should be used for annotation, for example:
   - `RNA_snn_res.0.9`
9. First draw the cluster UMAP and review it with the user.
10. Then compute cluster-level marker genes and marker isoforms before asking the user whether they want automatic annotation or manual annotation.

## Workflow

### 1. Load Libraries And Seurat Object

Within the persistent R session:

- `source("/work/.pi/scripts/annotation_helpers.r")`
- load the helper-managed libraries
- read the Seurat object RDS

### 2. Draw Cluster UMAP

Before annotation, first draw the cluster UMAP so the user can inspect the cluster structure.

Prefer using:

- `plot_cluster_umap()`

This helper uses the same `SCP::CellDimPlot` style as the QC workflow.

### 3. Compute Cluster-Level Marker Features

Compute:

- cluster marker genes from the `RNA` assay
- cluster marker isoforms from the isoform assay

Recommended project thresholds from the source text:

- `logfc.threshold = 0.5`
- `min.pct = 0.20`
- `only.pos = TRUE`
- keep markers with `p_val_adj < 0.05`

Save:

- cluster marker genes table
- cluster marker isoforms table

Prefer using:

- `compute_cluster_markers()`

### 4. Ask The User To Choose Annotation Mode

After cluster markers are available, ask the user whether they want:

- automatic annotation
- manual annotation

Do not assume one mode by default.

### 5. Automatic Annotation

If the user chooses automatic annotation, use the ScType-based workflow:

- load `gene_sets_prepare.R`
- load `sctype_score_.R`
- let the user specify:
  - database
  - tissue
  - optional gene-set removals
  - clustering metadata column
  - reduction name
  - metadata column prefix for the annotation result

Then:

- compute ScType scores
- assign one top cell type per cluster
- set low-confidence clusters to `Unknown`
- add the predicted labels into Seurat metadata
- draw the annotation UMAP

Prefer using:

- `perform_sctype_annotation()`

### 6. Manual Annotation

If the user chooses manual annotation, annotate cluster by cluster in the same persistent R session.

Requirements:

- ask the user to provide the cell type label for each cluster
- rename cluster identities accordingly
- save the final labels into a metadata column such as `celltype`
- draw the celltype UMAP after annotation

Prefer using:

- `rename_clusters_manually()`
- `plot_celltype_umap()`

### 7. Feature Plot Support During Annotation

During either automatic or manual annotation, if the user wants to inspect specific genes or isoforms, use:

- `../draw-features/SKILL.md`

Use that skill inside the same annotation session context to support annotation decisions.

### 8. Celltype-Level Differential Analysis

After annotation is finalized, perform celltype-level differential analysis for:

- genes
- isoforms

Using the project pattern in the source text:

- `RunDEtest`
- `VolcanoPlot`

Save:

- celltype DE gene table
- celltype DE isoform table

Prefer using:

- `run_celltype_de()`
- `plot_volcano()`

### 9. Heatmaps

Generate heatmaps for:

- celltype differential genes
- celltype differential isoforms

The source text uses:

- `AnnotateFeatures`
- `FeatureHeatmap`

For genes, include feature annotation and optional pathway annotations where appropriate.

Prefer using:

- `plot_celltype_heatmap()`

### 10. Enrichment Analysis

Run celltype-level enrichment analysis for differential genes using the project pattern:

- `RunEnrichment`
- `EnrichmentPlot`

Supported plot styles in the source text include:

- bar
- wordcloud
- feature wordcloud
- network

Prefer using:

- `run_celltype_enrichment()`
- `plot_enrichment()`

### 11. Cell Composition Summaries

Generate cell composition summaries such as:

- celltype by sample
- sample by celltype

The source text uses:

- `CellStatPlot`

Prefer using:

- `plot_cell_statistics()`

## Mandatory Rules

1. Keep the whole workflow in one persistent `tmux` R session until the user says to exit.
2. At the beginning of the session, first `source("/work/.pi/scripts/annotation_helpers.r")`.
3. Prefer calling helper functions from `annotation_helpers.r` instead of rewriting the same blocks manually in the terminal.
4. First show the cluster UMAP before annotation decisions.
5. First compute cluster-level marker genes and isoforms before asking the user to choose annotation mode.
6. Do not silently choose automatic or manual annotation for the user.
7. During manual annotation, do not rename clusters without explicit user-provided labels.
8. During annotation, if the user asks for gene or isoform feature plots, use the `draw-features` skill rather than inventing a separate plotting workflow.
9. After annotation is finalized, save the annotation labels back into Seurat metadata before running celltype-level DE analysis.
10. Before any substantial command or code block is executed, show the exact intended code or command and wait for confirmation.
11. Do not claim annotation is complete until the expected label UMAPs, DE tables, and summary plots exist and are non-empty.

## Expected Inputs

Typical inputs for this skill may include:

- merged Seurat-object RDS
- output directory
- clustering metadata column
- optional ScType database and tissue
- optional feature lists to inspect during annotation

## Expected Deliverables

For a normal successful annotation workflow, produce only the artifacts relevant to the chosen path, for example:

- cluster UMAP
- cluster marker gene table
- cluster marker isoform table
- automatic or manual annotation UMAP
- celltype metadata stored in the Seurat object
- celltype DE gene table
- celltype DE isoform table
- volcano plots
- heatmaps
- enrichment plots
- cell composition plots

## Response Style

Be operational and session-aware:

- state that annotation is being done in one persistent `tmux` R session
- say which Seurat object and clustering column are being used
- ask the user to choose automatic or manual annotation only after markers are ready
- if the user requests feature plots during annotation, keep the same overall annotation session alive
