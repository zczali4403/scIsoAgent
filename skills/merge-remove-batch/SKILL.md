---
name: merge-remove-batch
description: Merges per-sample Seurat objects after QC and doublet removal, generates PCA and Harmony elbow plots for user-confirmed dimension choices, previews batch structure in UMAP, and optionally performs Harmony-based batch correction before saving the final merged object.
---

# Merge Remove Batch

Use this skill when the user has finished per-sample QC and doublet removal and now wants to merge samples into one downstream Seurat object, with an explicit user decision on whether batch correction is needed.

This skill is for post-QC sample merging and optional batch correction only. It does not cover empty-droplet removal, per-sample QC, cell type annotation, marker discovery, or trajectory analysis.

## Use This Skill When

- the user has one QC-completed singlet Seurat object per sample
- the user wants to merge multiple samples into one Seurat object
- the user wants to inspect whether batch effects are visible before deciding on integration
- the user wants Harmony-based batch correction after merging

## Project Source

Execution source of truth:

- `../../scripts/merge_remove_batch.r`

Prefer calling that script directly.

Related upstream skill:

- `../qc-pipeline/SKILL.md`

## Current Scope

This skill currently includes:

1. collecting one singlet Seurat-object RDS per sample after QC
2. merging all per-sample objects into one Seurat object
3. running preprocessing to PCA
4. generating a PCA elbow plot
5. asking the user to confirm `npc`
6. generating a no-batch-correction preview UMAP by sample and by cluster
7. asking the user whether batch correction is needed
8. if Harmony is requested, running Harmony and generating a Harmony elbow plot
9. asking the user to confirm `harmony_npc`
10. running final neighbors, clusters, and UMAP
11. saving the final merged Seurat object and final UMAP plots

## Mandatory Rules

1. Treat this skill as a staged workflow, not as a one-shot merge command.
2. Run this R-based skill in the `downstream` conda environment.
3. First resolve the per-sample singlet Seurat-object RDS inputs.
4. The per-sample RDS paths may be provided by the user directly, or auto-detected from prior QC outputs and then confirmed by the user.
5. Do not silently assume auto-detected inputs without showing the resolved sample-to-file mapping and asking the user to confirm it.
6. Prefer calling `../../scripts/merge_remove_batch.r` with an explicit `--stage`.
7. When using this script for gene-level merging, pass `--object_label=gene`.
8. Before the preview UMAP stage, first generate the PCA elbow plot and ask the user to confirm `npc`.
9. After the preview UMAP is generated, ask the user whether batch correction is needed. Do not run Harmony by default.
10. If Harmony is requested, first generate the Harmony elbow plot and ask the user to confirm `harmony_npc` before `FindNeighbors`.
11. When Harmony is used, use the script's Harmony workflow as the execution source of truth rather than composing ad hoc integration code.
12. Before executing any substantial command or script block, show the exact final command and wait for confirmation.
13. Do not claim completion until the expected RDS, UMAP PDFs, and summary tables exist and are non-empty.

## Entry Procedure

At the start of a merge request:

1. Resolve the per-sample QC-completed singlet Seurat-object RDS files:
   - either ask the user to provide them directly
   - or auto-detect likely QC outputs and ask the user to confirm the resolved paths
2. Derive the merge output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
3. Confirm the output prefix for the merged object.
4. Use `09_merge_remove_batch` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
5. Optionally infer sample names from the input file names, then ask the user to confirm them if they will be used as the sample labels.
6. Validate that all declared input files exist before drafting commands.

Only after those are clear should the pipeline enter the first merge stage.

## Workflow Stages

### 1. PCA Elbow Plot

Merge the per-sample objects and run:

- `NormalizeData`
- `FindVariableFeatures`
- `ScaleData`
- `RunPCA`

Then generate the PCA elbow plot.

Recommended script stage:

- `--stage=elbow`

### 2. User-Confirmed PCA Dimensions

After the PCA elbow plot is generated, ask the user to confirm:

- `npc`

Do not proceed to UMAP preview before this confirmation.

### 3. No-Batch-Correction Preview

Using the user-confirmed `npc`, run:

- `FindNeighbors(reduction = "pca", dims = 1:npc)`
- `FindClusters`
- `RunUMAP(reduction = "pca", dims = 1:npc)`

Then generate:

- UMAP colored by sample
- UMAP colored by cluster

Recommended script stage:

- `--stage=preview`

### 4. User Decision On Batch Correction

After the preview UMAP is generated, ask the user whether batch correction is needed.

If the answer is no:

- proceed directly to finalization with `--remove_batch=false`

If the answer is yes:

- proceed to Harmony elbow plotting

### 5. Harmony Elbow Plot

If batch correction is requested, run Harmony first using the PCA dimensions selected above, then generate the Harmony elbow plot.

Recommended script stage:

- `--stage=harmony-elbow`

After the Harmony elbow plot is generated, ask the user to confirm:

- `harmony_npc`

Do not proceed to Harmony-based `FindNeighbors` before this confirmation.

### 6. Finalization

If batch correction is not requested:

- save the merged Seurat object after the non-Harmony neighbors/clusters/UMAP workflow

If batch correction is requested:

- run `FindNeighbors`, `FindClusters`, and `RunUMAP` on the Harmony reduction using the user-confirmed `harmony_npc`
- save the batch-corrected merged Seurat object
- save final UMAPs colored by sample and by cluster

Recommended script stage:

- `--stage=finalize`

## Expected Inputs

Typical inputs for this skill may include:

- one QC-completed singlet Seurat-object RDS per sample
- output directory
- output prefix
- optional explicit sample names

## Script Interface

Preferred execution pattern uses staged commands.

PCA elbow plot:

```bash
conda run -n downstream Rscript ../../scripts/merge_remove_batch.r \
  --stage=elbow \
  --input_rds_files=/abs/path/to/sample1_singlets.rds,/abs/path/to/sample2_singlets.rds \
  --output_dir=/abs/path/to/merge_remove_batch \
  --output_prefix=merged_samples \
  --sample_names=sample1,sample2 \
  --object_label=gene
```

Preview UMAP before batch correction:

```bash
conda run -n downstream Rscript ../../scripts/merge_remove_batch.r \
  --stage=preview \
  --input_rds_files=/abs/path/to/sample1_singlets.rds,/abs/path/to/sample2_singlets.rds \
  --output_dir=/abs/path/to/merge_remove_batch \
  --output_prefix=merged_samples \
  --sample_names=sample1,sample2 \
  --object_label=gene \
  --npc=30 \
  --cluster_res=0.5
```

Harmony elbow plot:

```bash
conda run -n downstream Rscript ../../scripts/merge_remove_batch.r \
  --stage=harmony-elbow \
  --input_rds_files=/abs/path/to/sample1_singlets.rds,/abs/path/to/sample2_singlets.rds \
  --output_dir=/abs/path/to/merge_remove_batch \
  --output_prefix=merged_samples \
  --sample_names=sample1,sample2 \
  --object_label=gene \
  --npc=30 \
  --batch_method=harmony
```

Finalize without batch correction:

```bash
conda run -n downstream Rscript ../../scripts/merge_remove_batch.r \
  --stage=finalize \
  --input_rds_files=/abs/path/to/sample1_singlets.rds,/abs/path/to/sample2_singlets.rds \
  --output_dir=/abs/path/to/merge_remove_batch \
  --output_prefix=merged_samples \
  --sample_names=sample1,sample2 \
  --object_label=gene \
  --npc=30 \
  --cluster_res=0.5 \
  --remove_batch=false
```

Finalize with Harmony:

```bash
conda run -n downstream Rscript ../../scripts/merge_remove_batch.r \
  --stage=finalize \
  --input_rds_files=/abs/path/to/sample1_singlets.rds,/abs/path/to/sample2_singlets.rds \
  --output_dir=/abs/path/to/merge_remove_batch \
  --output_prefix=merged_samples \
  --sample_names=sample1,sample2 \
  --object_label=gene \
  --npc=30 \
  --harmony_npc=20 \
  --cluster_res=0.5 \
  --remove_batch=true \
  --batch_method=harmony
```

## Expected Deliverables

For a normal successful run, produce:

- `<output_prefix>_sample_summary.tsv`
- `<output_prefix>_run_summary.tsv`
- `<output_prefix>_elbowplot.pdf`
- `<output_prefix>_preview_by_sample_umap.pdf`
- `<output_prefix>_preview_by_cluster_umap.pdf`
- when Harmony is requested:
  - `<output_prefix>_harmony_elbowplot.pdf`
- final outputs:
  - `<output_prefix>.rds`
  - `<output_prefix>_by_sample_umap.pdf`
  - `<output_prefix>_by_cluster_umap.pdf`

## Response Style

Be operational and staged:

- say which per-sample RDS files are being merged
- say which `npc` and `cluster_res` are in use
- after preview, explicitly ask whether batch correction is needed
- if Harmony is chosen, explicitly ask for `harmony_npc` after showing the Harmony elbow plot
- show the exact command before execution
