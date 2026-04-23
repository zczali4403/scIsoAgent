---
name: exploring-isoforms-of-interest
description: Explores isoforms of interest in one persistent tmux R session, including pseudobulk-based isoform counting per gene, ranking genes by isoform complexity, feature-plot exploration of selected genes, celltype-level violin and heatmap summaries for selected isoforms, and IsoVis input export.
---

# Exploring Isoforms Of Interest

Use this skill when the user wants focused downstream exploration of selected genes or isoforms after a Seurat object with an isoform assay is already available.

This skill is for exploratory isoform analysis only. It does not cover QC, cell type annotation itself, trajectory analysis, or DTU testing.

## Mandatory Session Rule

All analysis in this skill must be performed in one persistent interactive `tmux` R session.

Do not split this workflow across repeated `Rscript` calls.
Do not exit the R session after one substep.
Keep the same `tmux` R terminal alive until the user explicitly says to exit this exploration workflow.

## Use This Skill When

- the user wants to quantify how many isoforms each gene has
- the user wants to identify genes with the largest number of transcript isoforms
- the user wants to inspect one gene's transcript isoforms in feature plots
- the user wants violin, dot, or heatmap summaries of selected isoforms across cell types
- the user wants to prepare pseudobulk-format files for IsoVis

## Project Source Pattern

This skill is primarily an interactive workflow skill executed inside one persistent R session.

It relies on a Seurat object that already contains:

- a gene assay
- an isoform assay
- cell identity or annotation metadata such as `celltype` or `sctype_db`

It may also use the gene-symbol mapping file created earlier by `resource-setup`, especially when gene labels in rank plots should be shown as gene symbols instead of ENSG IDs.

Likely libraries needed in the session include:

- `Seurat`
- `SCP`
- `tidyverse`
- `ggplot2`
- `ggrepel`
- `Matrix`

Primary local helper script:

- `../scripts/exploring_isoforms_helpers.r`

At the start of the exploration session, prefer:

```r
source("../scripts/exploring_isoforms_helpers.r")
```

Then call helper functions from that script instead of repeatedly rewriting the same exploration code in the terminal.

## Entry Procedure

At the start of an isoform exploration request:

1. Confirm which conda environment should be used.
2. Start one persistent `tmux` R session and keep it open for the whole workflow.
3. Resolve the Seurat-object RDS path to explore.
4. In that session, `source("../scripts/exploring_isoforms_helpers.r")`.
5. Read the Seurat object once in that session.
6. Confirm which metadata column should be used for grouping, for example:
   - `celltype`
   - `sctype_db`
7. Derive the output directory for exploration artifacts from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
8. Use `14_explore_isoforms` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
9. If the user wants human-readable gene labels in ranked-gene outputs, resolve the gene-symbol mapping file created by `resource-setup`.
10. Ask which part of the exploration workflow the user wants first:
   - global isoform-count overview
   - top genes by isoform complexity
   - one gene's isoform exploration
   - IsoVis export

## Workflow

### 1. Load Seurat Object

Within the persistent R session:

- `source("../scripts/exploring_isoforms_helpers.r")`
- read the Seurat RDS
- ensure the grouping metadata column is available
- if needed, copy current identities into a metadata column such as `celltype`

Prefer using:

- `read_isoform_exploration_object()`

### 2. Global Isoform Count Per Gene

Aggregate isoform expression across groups using the isoform assay, then:

- split transcript and gene identifiers
- count the number of distinct isoforms per gene
- bin genes into categories such as:
  - `1`
  - `2-3`
  - `4-5`
  - `6+`
- plot the percentage of genes in each category

This is the high-level summary of isoform complexity across genes.

Prefer using:

- `aggregate_isoform_expression()`
- `split_isoform_gene_ids()`
- `compute_isoform_count_per_gene()`
- `categorize_isoform_count_per_gene()`
- `summarize_isoform_categories()`
- `plot_isoform_count_summary()`

### 3. Genes With The Largest Number Of Isoforms

Using the same pseudobulk-derived table:

- count unique transcripts per gene
- rank genes by isoform count
- write the full ranked table to CSV
- highlight the top `n` genes in a rank plot

The user should be allowed to choose `n`.

If the user wants the highlighted top-gene labels shown as gene symbols rather than ENSG IDs, use the gene-symbol mapping file from `resource-setup`.

Prefer using:

- `read_gene_symbol_mapping()`
- `rank_genes_by_isoform_count()`
- `plot_top_isoform_rich_genes()`

### 4. Explore One Gene's Isoforms

Ask the user to provide a gene of interest.

Then:

- find matching isoform features for that gene in the isoform assay
- compute total expression per isoform
- rank isoforms by total expression
- show the ranked isoform list to the user

Do not assume the number of isoforms to plot is always 12. Let the user choose how many top isoforms to visualize.

Prefer using:

- `rank_isoforms_for_gene()`
- `write_ranked_isoforms_for_gene()`

### 5. Feature Plots For Selected Isoforms

For the selected gene's isoforms:

- draw feature plots for the chosen isoforms on the chosen reduction
- keep the session open after each batch so the user can continue exploring other genes or isoforms

Prefer using:

- `SCP::FeatureDimPlot(..., theme_use = "theme_blank")`

If detailed feature plotting control is needed, use:

- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/draw-features/SKILL.md`

For quick plotting directly in this exploration workflow, helper support is also available via:

- `plot_selected_isoforms_featureplot()`

### 6. Celltype-Level Violin, Dot, And Heatmap Summaries

For user-selected isoforms:

- draw celltype-level expression summaries using `FeatureStatPlot`
- draw grouped isoform heatmaps using `GroupHeatmap`

The user should choose which isoforms are included in these summaries.

Prefer using:

- `plot_selected_isoforms_stat()`
- `plot_selected_isoforms_heatmap()`

### 7. IsoVis Export

Prepare a pseudobulk-format table for IsoVis by:

- exporting the pseudobulk isoform expression table to CSV

Then explicitly tell the user:

- IsoVis also requires a GTF file in addition to this expression CSV

Prefer using:

- `export_isovis_pseudobulk()`

## Mandatory Rules

1. Keep the whole workflow in one persistent `tmux` R session until the user says to exit.
2. At the beginning of the session, first `source("../scripts/exploring_isoforms_helpers.r")`.
3. Prefer calling helper functions from `exploring_isoforms_helpers.r` instead of rewriting the same blocks manually in the terminal.
4. Do not reload the Seurat object for every subtask.
5. Do not assume the grouping column; confirm it.
6. Do not assume one fixed gene of interest; ask the user.
7. Do not assume one fixed number of isoforms to visualize; ask the user how many top isoforms they want to inspect.
8. If the user wants gene labels shown as gene symbols in isoform-rich gene rank plots, prefer using the mapping file generated by `resource-setup`.
9. When that mapping file is available, do not leave top-gene labels as ENSG IDs by default.
10. If the user wants additional genes or isoforms after one round of plots, continue in the same `tmux` R session.
11. If the user wants feature plots, use the `draw-features` skill when its more structured plotting workflow is helpful.
12. Before any substantial command or code block is executed, show the exact intended code or command and wait for confirmation.
13. Do not claim exploration is complete until the expected plots, tables, or IsoVis export files exist and are non-empty.

## Expected Inputs

Typical inputs for this skill may include:

- Seurat-object RDS with isoform assay
- conda environment name
- grouping metadata column
- output directory
- optional gene-symbol mapping file from `resource-setup`
- optional gene of interest
- optional top `n` for ranking or plotting

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- isoform-count-per-gene summary plot
- ranked genes by isoform count CSV
- top-gene isoform complexity rank plot
- one gene's ranked isoform expression table
- feature plots for selected isoforms
- violin or dot summaries across cell types
- grouped heatmaps for selected isoforms
- IsoVis pseudobulk export CSV

## Response Style

Be operational and session-aware:

- state that exploration is being done in one persistent `tmux` R session
- say which Seurat object and grouping column are being used
- ask the user which exploration branch they want first
- keep the session alive after one gene or one plot batch so the user can continue exploring
