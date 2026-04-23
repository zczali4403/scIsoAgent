---
name: remove-droplets
description: Runs project-specific empty-droplet removal for long-read single-cell RNA-seq count matrices using the remove_droplets.r script, producing QC plots, summary tables, and Seurat objects before and after droplet filtering.
---

# Remove Droplets

Use this skill when the downstream workflow needs empty-droplet detection and removal before later QC, ambient RNA correction, clustering, or annotation.

This skill is for empty-droplet removal only. It does not cover ambient RNA correction, doublet detection, clustering, cell type annotation, or trajectory analysis.
When multiple per-sample matrices are present, run this step once per sample rather than once on the merged multi-sample matrix.

## Use This Skill When

- the user wants to remove empty droplets
- the user wants to run `emptyDrops`
- the user has a gene count matrix ready for empty-droplet calling
- the user has one per-sample gene matrix, or a set of per-sample matrices produced by `split-gene-count-by-sample`
- the downstream workflow needs filtered Seurat objects after droplet cleanup

## Project Source

Execution source of truth:

- `../scripts/remove_droplets.r`

Prefer calling that script directly.

## Core Outputs

This skill currently produces:

1. a PDF report with:
   - barcode rank plot
   - UMAP colored by droplet classification
   - UMAPs for `nCount_RNA` and `nFeature_RNA`
   - a summary table
2. a CSV summary file
3. a Seurat object before droplet removal
4. a Seurat object after removing empty droplets

## Mandatory Rules

1. First confirm all required input files before drafting commands.
2. Prefer calling the project script directly:
   - `../scripts/remove_droplets.r`
3. Collect the script parameters explicitly before execution.
4. Before execution, first ask the user which conda environment should be used for this script.
5. Do not assume the current shell environment is the correct R environment.
6. Use the main gene count matrix as the direct input to `emptyDrops` in this project script.
7. Always ask the user to confirm the `fdr_threshold` and `lower` cutoff before execution, while providing the project default values as the recommended starting point.
8. Prefer detecting input matrices automatically from the immediately preceding `split-gene-count-by-sample` outputs when those outputs exist in the current project context.
9. If automatic detection is missing, ambiguous, or clearly points to the wrong files, ask the user to provide the paths manually.
10. When auto-detecting split outputs, first look for per-sample files matching:
    - `*_geneSymbol_gene_count.csv`
11. Treat those files as the preferred input candidates for `gene_count_file`.
12. If multiple per-sample matrices are detected, iterate over them and run `remove_droplets.r` once per sample.
13. Use a sample-specific `output_seurat_file` and output directory or subdirectory for each run so outputs do not overwrite each other.
14. Before execution, show the exact final command or command set and wait for confirmation.
15. Do not claim completion until the expected PDF, CSV, and RDS outputs exist and are non-empty for every sample that was processed.
16. Reuse the already confirmed workflow-level results root and derive the remove-droplets output subdirectory automatically; do not ask the user to type a fresh output directory unless they explicitly want to override it.
17. Use `07_remove_droplets` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Entry Procedure

At the start of a droplet-removal request:

1. Confirm the main gene count matrix path:
   - first try to detect it from the most recent `split-gene-count-by-sample` result
   - when detecting automatically, prefer files named `*_geneSymbol_gene_count.csv`
   - otherwise ask for `gene_count_file`
2. If multiple per-sample matrices are detected, report the full sample list and state that the script will be run once per sample.
3. Derive the output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path:
   - `output_path`
4. Confirm how sample-specific outputs should be named:
   - by default, use one `output_seurat_file` per sample derived from the sample name
5. Ask which conda environment should be used.
6. Ask the user to confirm the threshold settings, and provide these default values as the recommended starting point:
   - `fdr_threshold=0.001`
   - `lower=500`
7. Optionally collect a project label for Seurat object creation:
   - `project_name`
8. Validate that all declared inputs exist before building commands.

If automatic input detection is used, explicitly report which files were resolved before showing the final command.
If multiple samples will be processed, explicitly report the per-sample command pattern before execution.

## Script Interface

Preferred execution pattern:

```bash
conda run -n <env_name> Rscript ../scripts/remove_droplets.r \
  --gene_count_file=/abs/path/to/gene_count.csv \
  --output_path=/abs/path/to/output_dir \
  --output_seurat_file=sample_name \
  --fdr_threshold=0.001 \
  --lower=500 \
  --project_name=remove_droplets
```

Required parameters:

- `env_name`
- `gene_count_file`
- `output_path`
- `output_seurat_file`

Optional parameters:

- `fdr_threshold`
- `lower`
- `project_name`

## Expected Deliverables

For a normal successful run, produce:

- one set of outputs per sample:
  - `<output_seurat_file>_plots.pdf`
  - `<output_seurat_file>_summary.csv`
  - `with_empty_<output_seurat_file>.rds`
  - `removed_empty_<output_seurat_file>.rds`

## Response Style

Be operational and explicit:

- say which inputs are being used
- say which thresholds are being applied
- say which output files are expected next
- show the exact command before execution
- do not mix empty-droplet removal with later downstream QC modules inside this skill
