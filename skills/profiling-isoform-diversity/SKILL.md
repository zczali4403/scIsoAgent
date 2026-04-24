---
name: profiling-isoform-diversity
description: Profiles isoform diversity in one persistent tmux R session using a Seurat object with an isoform assay, including celltype-specific subsetting, entropy-based diversity analysis with find_diversity.R, entropy heatmaps, high- and low-entropy gene summaries, gene-specific entropy plots, and per-cell isoform proportion pie charts.
---

# Profiling Isoform Diversity

Use this skill when the user wants downstream isoform diversity analysis after a Seurat object with an isoform assay is already available.

This skill is for entropy-based isoform diversity profiling only. It does not cover QC, cell type annotation itself, trajectory analysis, or DTU.

## Mandatory Session Rule

All analysis in this skill must be performed in one persistent interactive `tmux` R session.

Do not split this workflow across repeated `Rscript` calls.
Do not exit the R session after one substep.
Keep the same `tmux` R terminal alive until the user explicitly says to exit this diversity workflow.

## Project Source Pattern

This skill depends on these project scripts:

- `../../scripts/find_diversity.R`
- `../../scripts/profiling_isoform_diversity_helpers.r`

At the start of the diversity session, prefer:

```r
source("../../scripts/find_diversity.R")
source("../../scripts/profiling_isoform_diversity_helpers.r")
```

## Use This Skill When

- the user wants entropy-based isoform diversity analysis
- the user wants to study isoform diversity in one selected cell type
- the user wants entropy heatmaps across genes and cells
- the user wants to identify genes with the highest or lowest isoform entropy
- the user wants to inspect selected genes in more detail
- the user wants per-cell isoform proportion pie charts for a chosen gene

## Entry Procedure

At the start of an isoform-diversity request:

1. Confirm which conda environment should be used.
2. Start one persistent `tmux` R session and keep it open for the whole workflow.
3. Resolve the Seurat-object RDS path to analyze.
4. In that session, first `source("../../scripts/find_diversity.R")`, then `source("../../scripts/profiling_isoform_diversity_helpers.r")`.
5. Read the Seurat object once in that session.
6. Confirm which metadata column should be used for celltype selection, for example:
   - `celltype`
   - `sctype_db`
7. Ask the user which cell type should be analyzed first.
8. Derive the output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
9. Use `14_isoform_diversity` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
10. Resolve gene-symbol mapping configuration for entropy row labels:
   - do not hardcode a fixed absolute dictionary path in commands
   - prefer user-provided dictionary path if given
   - otherwise auto-discover under `<output_root>/resource_setup/*isoform_gene_symbol_dict.csv`
   - if no dictionary is found, continue with `gene_id` labels and explicitly report that to the user

## Workflow

### 1. Load Seurat Object And Subset One Cell Type

Within the persistent R session:

- read the Seurat-object RDS, preferably with `read_isoform_diversity_object()`
- ask the user which cell type should be analyzed
- subset the Seurat object to that cell type, preferably with `subset_diversity_celltype()`

The source text uses a pattern like:

```r
RG_seu_obj <- subset(seu_obj, subset = sctype_db == "Radial glial cells")
```

Do not assume the target cell type; always ask the user.

### 2. Run Isoform Diversity Analysis

On the celltype-specific subset:

- set the default assay as needed
- define the genes to analyze
- prefer `run_isoform_diversity()` from `profiling_isoform_diversity_helpers.r`, which wraps `find_diversity()`
- when preparing entropy row labels, prefer gene-symbol mapping when a dictionary is available; otherwise keep gene_id labels

The source text uses:

- `min_counts_per_cell = 10`
- `assay = "iso"`
- `genes = rownames(RG_seu_obj)`

Save the resulting diversity object, for example:

- `res.rds`

### 3. Entropy Heatmap

Using `res$normalized_entropy`:

- prefer `plot_entropy_heatmap()`
- draw a heatmap where rows are genes and columns are cells
- use entropy as the color scale
- hide row and column names by default

The source text uses `pheatmap` and suggests colors such as:

- `"#3288BD"`
- `"white"`
- `"#D53E4F"`

### 4. High- and Low-Entropy Genes

Compute median normalized entropy per gene, then:

- prefer `summarize_entropy_by_gene()`
- identify genes with the lowest entropy
- identify genes with the highest entropy
- draw summary plots for low-entropy and high-entropy genes, preferably with `plot_high_low_entropy_genes()`
- display the global median as a reference line

The source text explicitly selects:

- low-entropy genes from the smallest 25 nonzero median-entropy genes

### 5. Gene-Specific Entropy Inspection

Ask the user which genes should be inspected further.

Then:

- subset the normalized entropy matrix to those genes, preferably with `prepare_gene_entropy_long()`
- reshape to long format
- map cells back to their celltype labels
- draw per-gene entropy plots using quasirandom points and median bars, preferably with `plot_gene_entropy()` or `plot_gene_entropy_batch()`

Do not assume one fixed set of genes; always ask the user.

### 6. Isoform Proportion Pie Charts For Selected Cells

For a user-selected gene:

- show isoform-expression proportion pie charts across a chosen number of cells
- ask the user to confirm:
  - the gene
  - `top_n`
  - `n_cells`
  - plotting layout such as `plot_col`

The source text requires a helper function:

- `plotIsoformPieCells()`

This function should draw, for each selected cell:

- the top `n` isoforms
- one aggregated `other` category

Use the provided `color_palette_21` for transcript colors.

## Mandatory Rules

1. Keep the whole workflow in one persistent `tmux` R session until the user says to exit.
2. At the beginning of the session, first `source("../../scripts/find_diversity.R")`, then `source("../../scripts/profiling_isoform_diversity_helpers.r")`.
3. Prefer helper functions from `profiling_isoform_diversity_helpers.r` for loading objects, subsetting cell types, running diversity analysis, entropy summaries, and pie-chart plotting instead of rewriting the same code inline.
4. Do not reload the Seurat object for every subtask.
5. Do not assume the target cell type; ask the user.
6. Do not assume the genes for detailed entropy inspection; ask the user.
7. Do not assume the settings for pie-chart visualization; ask the user.
8. Before any substantial command or code block is executed, show the exact intended code or command and wait for confirmation.
9. Do not claim diversity profiling is complete until the expected result RDS and requested plots exist and are non-empty.
10. Do not hardcode one fixed absolute gene-symbol dictionary path in this workflow.
11. Gene-symbol dictionary resolution priority must be:
    - user-provided dictionary path
    - auto-discovery under `<output_root>/resource_setup/*isoform_gene_symbol_dict.csv`
    - fallback to gene_id labels with an explicit user-facing note.

## Expected Inputs

Typical inputs for this skill may include:

- Seurat-object RDS with isoform assay
- conda environment name
- metadata column for celltype selection
- target cell type
- output directory
- optional selected genes
- optional pie-chart settings

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- diversity result RDS
- entropy heatmap
- high-entropy gene summary plot
- low-entropy gene summary plot
- selected-gene entropy plots
- per-cell isoform proportion pie charts

## Response Style

Be operational and session-aware:

- state that diversity analysis is being done in one persistent `tmux` R session
- say which Seurat object, metadata column, and target cell type are being used
- ask the user which genes they want to inspect further
- keep the session alive after one plot batch so the user can continue exploring
