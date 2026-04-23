---
name: qc-pipeline
description: Orchestrates project-specific downstream quality control for long-read single-cell RNA-seq after empty-droplet removal, including metric plotting, user-confirmed filtering thresholds, dimensional reduction, clustering, doublet detection, and QC summary outputs.
---

# QC Pipeline

Use this skill when the user wants downstream quality control after empty-droplet removal and before later annotation, marker discovery, or trajectory analysis.

This skill is the QC workflow entry point. It covers staged QC plotting, threshold confirmation, normalization, PCA, clustering, doublet detection, and QC summary reporting. It should orchestrate the sequence of user decisions rather than ask for every parameter up front.

## Use This Skill When

- the user wants downstream QC after empty-droplet removal
- the user wants QC metric plots before filtering
- the user wants gene-level filtering based on feature counts, read counts, or mitochondrial percentage
- the user wants elbow-plot-guided PCA dimension choice
- the user wants doublet detection after clustering

## Project Source

- `../scripts/QC_pipeline.r`

Use the R script as the execution source of truth.

Related upstream/downstream skill:

- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/remove-droplets/SKILL.md`
- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/merge-remove-batch/SKILL.md`

## Current Scope

This QC pipeline currently includes:

1. loading the Seurat object after empty-droplet removal
2. generating pre-filter QC plots
3. asking the user to confirm filtering thresholds
4. filtering cells
5. generating post-filter QC plots
6. normalization, variable feature selection, scaling, and PCA
7. generating the elbow plot
8. asking the user to confirm `npc`
9. neighbor graph construction
10. asking the user to confirm clustering resolution
11. clustering and UMAP
12. doublet parameter sweep and `pK` selection
13. asking the user to confirm `doublet_rate`
14. doublet detection and singlet retention
15. final QC plots and summary tables

## Mandatory Rules

1. Treat this skill as a staged QC workflow, not as a single one-shot command.
2. Do not ask for all later-stage parameters at the beginning.
3. First confirm which conda environment should be used.
4. The input Seurat-object RDS path may be provided by the user directly, or auto-detected from the most recent `remove-droplets` result and then confirmed by the user.
5. Do not silently assume an auto-detected RDS path without showing the resolved path and asking the user to confirm it.
6. First confirm the input Seurat object after empty-droplet removal.
7. Prefer calling `../scripts/QC_pipeline.r` with an explicit `--stage` rather than composing ad hoc inline R code.
8. Use the same initial empty-droplet-filtered Seurat-object RDS as the input to every QC stage.
9. Do not require the user to pass a stage-specific intermediate RDS from the previous step.
10. Before filtering, generate the pre-filter QC plots and then ask the user to confirm filtering thresholds.
11. After PCA, generate the elbow plot and then ask the user to confirm `npc`.
12. After neighbor construction, ask the user to confirm clustering resolution before clustering and UMAP.
13. After `pK` has been estimated from the parameter sweep, ask the user to confirm `doublet_rate` before running `doubletFinder`.
14. Before executing any substantial command or script block, show the exact final command and wait for confirmation.
15. Do not claim QC is complete until the expected plots, summary tables, and filtered Seurat outputs exist and are non-empty.
16. Reuse the already confirmed workflow-level results root and derive the QC output subdirectory automatically; do not ask the user to type a fresh output directory unless they explicitly want to override it.
17. Use `08_qc_pipeline` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Entry Procedure

At the start of a QC request:

1. Confirm the conda environment to use.
2. Resolve the Seurat object RDS path produced after empty-droplet removal:
   - either ask the user to provide it directly
   - or auto-detect a likely `removed_empty_*.rds` output from the most recent `remove-droplets` result and ask the user to confirm that resolved path
3. Derive the QC output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
4. Confirm the sample identifier or output prefix.
5. Validate that the input Seurat object exists before proposing commands.

Only after those are clear should the pipeline enter the first QC stage.

## Workflow Stages

### 1. Pre-Filter QC Plotting

Load the Seurat object after empty-droplet removal and generate:

- scatter plot of `nCount_RNA` vs `nFeature_RNA`
- violin plots for:
  - `nFeature_RNA`
  - `nCount_RNA`
  - `percent.mt`

This stage is for inspection before filtering. Do not apply thresholds yet.

Recommended script stage:

- `--stage=pre-filter`

### 2. User-Confirmed Filtering Thresholds

After the pre-filter plots are generated, ask the user to confirm:

- `max.features`
- `min.features`
- `min.counts`
- `max.counts`
- `MT`

Provide context from the generated plots before asking for these values.

### 3. Cell Filtering

Apply the confirmed thresholds to subset the Seurat object.

Then generate the post-filter violin plot for:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mt`

Recommended script stage:

- `--stage=filter-pca`

### 4. Standard Preprocessing

Run:

- `NormalizeData`
- `FindVariableFeatures`
- `ScaleData`
- `RunPCA`

Then generate the elbow plot.

### 5. User-Confirmed PCA Dimensions

After the elbow plot is generated, ask the user to confirm:

- `npc`

Do not proceed to neighbors or clustering before this confirmation.

### 6. Neighbor Graph, Clustering, and UMAP

Run:

- `FindNeighbors(dims = 1:npc)`

Then ask the user to confirm:

- `cluster_res`

Only after that should you run:

- `FindClusters`
- `RunUMAP`

Recommended script stage:

- `--stage=cluster-umap`

### 7. Doublet Detection

Run:

- `paramSweep`
- `summarizeSweep`
- `find.pK`

Use the resulting `BCmetric` maximum to derive `pK`.

Then ask the user to confirm:

- `doublet_rate`

Recommended default:

- `doublet_rate = 0.05`

Only after confirmation should you run the final `doubletFinder` stage and subset to singlets.

Recommended script stage:

- `--stage=doublet`

### 8. Final QC Outputs

Generate:

- cluster UMAP
- `nCount_RNA` UMAP
- `nFeature_RNA` UMAP
- doublet-classification UMAP
- doublet statistics table
- overall QC summary table

## Expected Inputs

Typical inputs for this skill may include:

- empty-droplet-filtered Seurat object (`removed_empty_*.rds`) path
  - either manually provided by the user
  - or auto-detected from the most recent `remove-droplets` output and then user-confirmed
- conda environment name
- QC output directory
- sample identifier

## Downstream Handoff

After this skill has been completed for each sample and singlets have been retained, the next multi-sample downstream step is typically:

- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/merge-remove-batch/SKILL.md`

Use that skill when the user wants to:

- merge all per-sample QC-completed singlet Seurat objects
- inspect sample mixing in a merged UMAP
- decide whether Harmony-based batch correction is needed

## Script Interface

Preferred execution pattern uses staged commands.

Pre-filter plotting:

```bash
conda run -n <env_name> Rscript ../scripts/QC_pipeline.r \
  --stage=pre-filter \
  --input_rds=/abs/path/to/removed_empty_sample.rds \
  --output_dir=/abs/path/to/QC_pipeline \
  --sample_id=sample_name
```

Filtering plus PCA:

```bash
conda run -n <env_name> Rscript ../scripts/QC_pipeline.r \
  --stage=filter-pca \
  --input_rds=/abs/path/to/removed_empty_sample.rds \
  --output_dir=/abs/path/to/QC_pipeline \
  --sample_id=sample_name \
  --min_features=200 \
  --max_features=8000 \
  --min_counts=500 \
  --max_counts=50000 \
  --mt_cutoff=20
```

Clustering and UMAP:

```bash
conda run -n <env_name> Rscript ../scripts/QC_pipeline.r \
  --stage=cluster-umap \
  --input_rds=/abs/path/to/removed_empty_sample.rds \
  --output_dir=/abs/path/to/QC_pipeline \
  --sample_id=sample_name \
  --min_features=200 \
  --max_features=8000 \
  --min_counts=500 \
  --max_counts=50000 \
  --mt_cutoff=20 \
  --npc=20 \
  --cluster_res=0.5
```

Doublet detection:

```bash
conda run -n <env_name> Rscript ../scripts/QC_pipeline.r \
  --stage=doublet \
  --input_rds=/abs/path/to/removed_empty_sample.rds \
  --output_dir=/abs/path/to/QC_pipeline \
  --sample_id=sample_name \
  --min_features=200 \
  --max_features=8000 \
  --min_counts=500 \
  --max_counts=50000 \
  --mt_cutoff=20 \
  --npc=20 \
  --cluster_res=0.5 \
  --doublet_rate=0.05
```

## Expected Deliverables

For a normal successful run, produce:

- `nCount_nFeature_corr.pdf`
- `QC_metrics_violin.pdf`
- `QC_metrics_after_qc_violin.pdf`
- `Elbowplot.pdf`
- `seurat_clusters_umap.pdf`
- `nCount_RNA_umap.pdf`
- `nFeature_RNA_umap.pdf`
- `DF_classification_umap.pdf`
- `stats_Doublets.pdf`
- `stats_summary.pdf`
- filtered Seurat object with singlets retained

## Interaction Rules

Use progressive confirmation:

- first show pre-filter plots, then ask for filter thresholds
- then show elbow plot, then ask for `npc`
- then ask for `cluster_res`
- then ask for `doublet_rate`

Do not front-load these into one large parameter questionnaire.

## Response Style

Be operational and explicit:

- say which QC stage you are in
- say which input object is being used
- say which plots or outputs are expected next
- show exact commands before execution
- keep the workflow staged around the user confirmation points defined above
