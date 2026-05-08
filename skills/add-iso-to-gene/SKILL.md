---
name: add-iso-to-gene
description: Adds the merged isoform count assay and isoform UMAP from a merged isoform Seurat object into a merged gene Seurat object, then saves a new combined Seurat RDS for joint downstream gene-plus-isoform analysis.
---

# Add Iso To Gene

Use this skill when the user already has one merged gene-level Seurat object and one merged isoform-level Seurat object, and wants to combine them into a single Seurat object for joint downstream analysis.

This skill is for combining existing merged gene and isoform Seurat objects only. It does not cover per-sample isoform object creation, isoform-object merging, gene-level QC, or annotation.

## Use This Skill When

- the user has a merged gene Seurat object
- the user has a merged isoform Seurat object
- the user wants to attach isoform counts as a new assay inside the gene object
- the user wants to preserve the isoform UMAP inside the combined object
- the user wants one final RDS for joint gene-plus-isoform downstream analysis

## Project Source

Execution source of truth:

- `/work/.pi/scripts/add_iso_to_gene.r`

Prefer calling that script directly.

Related upstream skills:

- `../merge-remove-batch/SKILL.md`
- `../merge-iso-obj/SKILL.md`

## Current Scope

This skill currently includes:

1. loading one merged gene Seurat RDS
2. loading one merged isoform Seurat RDS
3. finding overlapping cell barcodes
4. subsetting both objects to the shared cell set
5. adding the isoform count matrix to the gene object as a new assay
6. adding the isoform UMAP reduction to the gene object under a separate reduction name when available
7. saving a combined Seurat RDS
8. writing one summary table

## Mandatory Rules

1. Run this R-based skill in the `downstream` conda environment.
2. First resolve the merged gene Seurat RDS path.
3. First resolve the merged isoform Seurat RDS path.
4. Prefer calling `/work/.pi/scripts/add_iso_to_gene.r` directly rather than composing inline R code.
5. Use this step only after both the merged gene object and merged isoform object already exist.
6. Before execution, confirm the output RDS path.
7. If the user wants custom assay or reduction names, collect them before building the final command.
8. Before executing any substantial command or script block, show the exact final command and wait for confirmation.
9. Do not claim completion until the combined RDS and summary table exist and are non-empty.
10. Reuse the already confirmed workflow-level results root and place default outputs under `12_add_iso_to_gene` unless the user explicitly asks for a different subdirectory name.

## Entry Procedure

At the start of an add-isoform-to-gene request:

1. Resolve the merged gene object path:
   - `gene_rds`
2. Resolve the merged isoform object path:
   - `isoform_rds`
3. Confirm the combined output path:
   - `output_rds`
4. Optionally confirm custom names for:
   - `isoform_assay_name`
   - `isoform_reduction_name`
5. Validate that all declared inputs exist before drafting commands.

Only after those are clear should the script be run.

## Workflow

### 1. Load The Two Merged Objects

Read:

- one merged gene Seurat RDS
- one merged isoform Seurat RDS

### 2. Align Shared Cells

Find overlapping cell barcodes between the two merged objects and subset both objects to the shared cell set.

### 3. Add Isoform Assay

Take the isoform object's count matrix and add it into the gene object as a new assay.

Default assay name:

- `isoform`

### 4. Add Isoform UMAP

If the merged isoform object contains a `umap` reduction, add it into the gene object under a separate reduction name.

Default reduction name:

- `isoform_umap`

### 5. Save Outputs

Save:

- the combined Seurat RDS
- one summary table describing inputs and shared-cell counts

## Expected Inputs

Typical inputs for this skill may include:

- `gene_rds`
- `isoform_rds`
- `output_rds`
- optional `isoform_assay_name`
- optional `isoform_reduction_name`

## Script Interface

Preferred execution pattern:

```bash
conda run -n downstream Rscript /work/.pi/scripts/add_iso_to_gene.r \
  --gene_rds=/abs/path/to/merged_gene.rds \
  --isoform_rds=/abs/path/to/merged_isoform.rds \
  --output_rds=/abs/path/to/gene_with_isoform.rds \
  --isoform_assay_name=isoform \
  --isoform_reduction_name=isoform_umap
```

Required parameters:

- `gene_rds`
- `isoform_rds`
- `output_rds`

Optional parameters:

- `isoform_assay_name`
- `isoform_reduction_name`

## Expected Deliverables

For a normal successful run, produce:

- the combined Seurat object:
  - `<output_rds>`
- the run summary:
  - `<output_rds_basename>_summary.tsv`

## Response Style

Be operational and explicit:

- say which merged gene object is being used
- say which merged isoform object is being used
- say which assay and reduction names will be added
- show the exact command before execution
