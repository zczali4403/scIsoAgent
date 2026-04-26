---
name: generate-iso-obj
description: Builds per-sample isoform-level Seurat objects from transcript barcode-feature-matrix triplets, plots isoform QC metrics, subsets each isoform object to the cells retained in the matched gene-level QC Seurat object, and saves the per-sample isoform Seurat RDS outputs.
---

# Generate Iso Obj

Use this skill when the user wants to construct per-sample isoform-level Seurat objects after gene-level empty-droplet removal and QC have already determined which cells should be retained.

This skill is for creating and filtering isoform Seurat objects only. It does not cover gene-level QC, sample merging, cell type annotation, marker detection, or trajectory analysis.

## Use This Skill When

- the user has per-sample transcript or isoform barcode-feature-matrix files from upstream quantification
- the user already has matched per-sample gene-level Seurat RDS files after `remove-droplets` and `QC-pipeline`
- the user wants to keep only the barcodes retained by the gene-level QC objects
- the user wants per-sample isoform Seurat objects saved as RDS files

## Project Source

Execution source of truth:

- `../../scripts/generate_iso_obj.r`

Prefer calling that script directly.

Related upstream skills:

- `../qc-pipeline/SKILL.md`
- `../resource-setup/SKILL.md`

## Current Scope

This skill currently includes:

1. resolving one transcript barcode-feature-matrix triplet per sample
2. resolving one matched gene-level Seurat RDS per sample
3. inferring sample names from transcript barcode file names using the same convention as `resource_file_creation.r` unless the user explicitly provides them
4. creating one isoform Seurat object per sample
5. plotting `nFeature_RNA` and `nCount_RNA` violin plots for each sample before subsetting
6. subsetting each isoform object to the barcodes present in the matched gene-level Seurat object
7. saving one filtered isoform Seurat RDS per sample
8. writing one summary table across all processed samples

## Mandatory Rules

1. Run this R-based skill in the `downstream` conda environment.
2. First resolve the per-sample transcript barcode-feature-matrix inputs.
3. First resolve the matched per-sample gene-level Seurat RDS inputs that came after empty-droplet removal and QC.
4. The sample order across `barcodes_paths`, `features_paths`, `matrix_paths`, and `gene_rds_paths` must match exactly.
5. Prefer calling `../../scripts/generate_iso_obj.r` directly rather than composing inline R code.
6. When `sample_names` are not explicitly provided, infer them from `barcodes_paths` using the same rule as `../../scripts/resource_file_creation.r` so gene and isoform barcode prefixes stay consistent.
7. Do not silently assume auto-detected inputs without showing the resolved sample-to-file mapping and asking the user to confirm it.
8. Before execution, show the exact final command and wait for confirmation.
9. Do not claim completion until the expected violin plots, isoform RDS outputs, and summary table exist and are non-empty for every sample.

## Entry Procedure

At the start of an isoform-object generation request:

1. Resolve the transcript-level per-sample inputs:
   - `barcodes_paths`
   - `features_paths`
   - `matrix_paths`
2. Resolve the matched gene-level QC Seurat inputs:
   - `gene_rds_paths`
3. Derive the output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
4. Use `10_generate_iso_obj` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
5. If `sample_names` are not explicitly provided, infer them from `barcodes_paths`, then show the inferred sample names and ask the user to confirm them.
6. Validate that all declared inputs exist and that the sample list lengths match before proposing commands.

Only after those are clear should the script be run.

## Workflow

### 1. Create Per-Sample Isoform Seurat Objects

For each sample, read:

- the transcript barcode file
- the transcript feature file
- the transcript matrix file

Then:

- prepend the sample name to each barcode
- create a Seurat object from the transcript matrix

### 2. Plot Isoform QC Metrics

For each sample, generate a violin plot showing:

- `nFeature_RNA`
- `nCount_RNA`

These plots are generated before any barcode filtering against the gene-level Seurat object.

### 3. Subset By Gene-Level QC Cells

For each sample, read the matched gene-level Seurat RDS and keep only the overlapping cell barcodes in the isoform Seurat object.

This step ensures the final isoform object uses the same retained cell set as the gene-level downstream workflow.

### 4. Save Outputs

For each sample, save:

- the isoform violin plot PDF
- the subsetted isoform Seurat RDS

Also save one cross-sample summary table.

## Expected Inputs

Typical inputs for this skill may include:

- `barcodes_paths`
- `features_paths`
- `matrix_paths`
- `gene_rds_paths`
- `output_dir`
- optional `sample_names`
- optional `min_features`

## Script Interface

Preferred execution pattern:

```bash
conda run -n downstream Rscript ../../scripts/generate_iso_obj.r \
  --barcodes_paths=/abs/path/to/sample1.barcodes.tsv,/abs/path/to/sample2.barcodes.tsv \
  --features_paths=/abs/path/to/sample1.features.tsv,/abs/path/to/sample2.features.tsv \
  --matrix_paths=/abs/path/to/sample1.matrix.mtx,/abs/path/to/sample2.matrix.mtx \
  --gene_rds_paths=/abs/path/to/sample1_gene_singlets.rds,/abs/path/to/sample2_gene_singlets.rds \
  --output_dir=/abs/path/to/generate_iso_obj \
  --sample_names=sample1,sample2 \
  --min_features=0
```

Required parameters:

- `barcodes_paths`
- `features_paths`
- `matrix_paths`
- `gene_rds_paths`
- `output_dir`

Optional parameters:

- `sample_names`
- `min_features`

## Expected Deliverables

For a normal successful run, produce:

- one set of outputs per sample:
  - `<output_dir>/<sample_name>/<sample_name>_isoform_violinplot.pdf`
  - `<output_dir>/<sample_name>/<sample_name>_isoform_subset.rds`
- one cross-sample summary:
  - `<output_dir>/generate_iso_obj_summary.tsv`

## Response Style

Be operational and explicit:

- say which transcript triplets are being used for each sample
- say which gene-level RDS is being used to subset each isoform object
- show the inferred or user-provided sample names before execution
- show the exact command before execution
