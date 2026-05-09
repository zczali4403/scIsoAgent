---
name: isoform-classification
description: Runs SQANTI3 isoform classification from the terminal and then summarizes SQANTI structural categories, subcategories, coding status, and celltype-specific isoforms in one persistent tmux R session using the merged Seurat object and pseudobulk isoform summaries.
---

# Isoform Classification

Use this skill when the user wants downstream isoform classification after a Seurat object with an isoform assay is already available and a transcriptome GTF can be passed to SQANTI3.

This skill covers SQANTI3-based structural classification plus downstream celltype-level summary analysis. It does not cover QC, cell type annotation itself, trajectory analysis, or DTU.

## Mandatory Session Rule

This skill has two execution phases:

1. SQANTI3 classification is run as terminal commands.
2. All subsequent R analysis must continue in one persistent interactive `tmux` R session until the user explicitly says to exit.

Do not split the R analysis across repeated `Rscript` calls.
Do not exit the R session after one substep.

## Project Source Pattern

This skill uses:

- terminal commands for SQANTI3
- one persistent `tmux` R session for downstream summaries

Primary local helper script:

- `/work/.pi/scripts/isoform_classification_helpers.r`

At the start of the R analysis session, prefer:

```r
source("/work/.pi/scripts/isoform_classification_helpers.r")
```

Then call helper functions from that script instead of repeatedly rewriting the same classification-summary code in the terminal.

## Use This Skill When

- the user wants SQANTI3-based isoform classification
- the user wants to summarize isoform structural categories across cell types
- the user wants to compare structural subcategories across cell types
- the user wants to compare coding proportions across cell types
- the user wants to detect celltype-specific isoforms

## Entry Procedure

At the start of an isoform-classification request:

1. Use the `sqanti3` conda environment for the terminal SQANTI3 phase.
2. Use the `downstream` conda environment for the persistent R analysis phase.
3. Resolve the transcriptome GTF to classify with SQANTI3.
4. Resolve the reference annotation GTF.
5. Resolve the reference genome FASTA.
6. Confirm the merged Seurat-object RDS used for downstream R analysis.
7. Confirm which metadata column should be used for grouping, for example:
   - `celltype`
   - `sctype_db`
8. Derive the output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
9. Use `14_isoform_classification` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
10. Run the SQANTI3 preprocessing and classification commands first.
11. Then start one persistent `tmux` R session in the `downstream` environment and keep it open for all remaining analysis.
12. In that session, `source("/work/.pi/scripts/isoform_classification_helpers.r")`.

## Workflow

### 1. SQANTI3 Classification

Before entering the R session, first preprocess the GTF and run SQANTI3 from the terminal.

The source pattern is:

```bash
conda run -n sqanti3 bash -lc "awk '$7 == \"+\" || $7 == \"-\"' isoform_annotated.gtf > remove_unknownstrand.gtf"
conda run -n sqanti3 /opt/SQANTI3/sqanti3_qc.py remove_unknownstrand.gtf <reference_gtf> <genome_fasta>
```

Requirements:

- ask the user to provide:
  - the transcriptome GTF to classify
  - the reference annotation GTF
  - the reference genome FASTA
- remove entries with unknown strand before running SQANTI3
- use the resulting `classification.txt` file for downstream R analysis
- start the subsequent persistent R session in the `downstream` environment

### 2. Load SQANTI Output And Seurat Object

Inside the persistent `tmux` R session:

- read the SQANTI `classification.txt`
- read the Seurat-object RDS
- aggregate isoform expression across cell types from the isoform assay
- split transcript and gene identifiers from the isoform feature names

Prefer using:

- `read_sqanti_classification()`
- `read_isoform_classification_object()`
- `aggregate_isoform_pseudobulk()`
- `split_isoform_feature_ids()`
- `build_long_isoform_expression()`

### 3. Celltype-Level Gene and Isoform Counts

From the pseudobulk isoform table:

- count the number of distinct genes per cell type
- count the number of distinct isoforms per cell type
- draw the bar plot summarizing genes and isoforms per cell type

Prefer using:

- `summarize_celltype_gene_isoform_counts()`
- `plot_celltype_gene_isoform_counts()`

### 4. Structural Category, Subcategory, and Coding Summaries

Merge the pseudobulk isoform table with the SQANTI classification table, then summarize across cell types:

- `structural_category`
- `subcategory`
- `coding`

Generate the plots shown in the source text:

- isoforms by structural category across cell types
- isoforms by structural subcategory across cell types
- proportion of coding isoforms across cell types

Prefer using:

- `merge_sqanti_with_pseudobulk()`
- `summarize_sqanti_categories_by_celltype()`
- `plot_structural_category_summary()`
- `plot_structural_subcategory_summary()`
- `plot_coding_proportion_summary()`
- `plot_sqanti_classification_overview()`

### 5. Celltype-Specific Isoforms

Detect celltype-specific isoforms using user-confirmed thresholds.

The source text uses:

- `expression_cutoff`
- `max_expression_in_all_other_cells_types`

You must ask the user to confirm these thresholds before running this step.

Interpretation rule from the source text:

- an isoform is considered celltype-specific if it is expressed above the chosen cutoff in at most one cell type and remains below the chosen background threshold in all other cell types

Then:

- count the number of celltype-specific isoforms per cell type
- draw the bar plot summarizing unique isoforms per cell type

Prefer using:

- `find_celltype_specific_isoforms()`
- `count_celltype_specific_isoforms()`
- `plot_celltype_specific_isoform_counts()`

## Mandatory Rules

1. Run SQANTI3 first, then do the R-based analysis in one persistent `tmux` R session.
2. At the beginning of the R session, first `source("/work/.pi/scripts/isoform_classification_helpers.r")`.
3. Prefer calling helper functions from `isoform_classification_helpers.r` instead of rewriting the same blocks manually in the terminal.
4. Do not assume the Seurat grouping column; confirm it.
5. Do not assume the thresholds for celltype-specific isoforms; ask the user.
6. Before any substantial command or code block is executed, show the exact intended command or code and wait for confirmation.
7. Do not claim classification is complete until the SQANTI output and the expected summary plots exist and are non-empty.

## Expected Inputs

Typical inputs for this skill may include:

- transcriptome GTF for SQANTI3
- reference annotation GTF
- reference genome FASTA
- SQANTI `classification.txt`
- merged Seurat-object RDS with isoform assay
- grouping metadata column
- output directory

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- strand-filtered GTF used for SQANTI3
- SQANTI `classification.txt`
- gene and isoform count summary across cell types
- structural category summary plots
- structural subcategory summary plots
- coding proportion plot
- celltype-specific isoform count plot

## Response Style

Be operational and session-aware:

- state when the workflow is in the terminal SQANTI step versus the persistent R analysis step
- say which Seurat object and grouping column are being used
- ask for thresholds before the celltype-specific isoform step
- keep the same `tmux` R session alive until the user explicitly says to exit
