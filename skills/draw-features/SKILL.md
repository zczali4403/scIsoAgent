---
name: draw-features
description: Draws gene or isoform feature expression on a chosen reduction from a Seurat RDS using SCP::FeatureDimPlot, and for repeated plotting must keep one R session open so the object is loaded once and the user can iteratively request more plots.
---

# Draw Features

Use this skill when the user wants to visualize one or more genes or isoforms on a Seurat reduction such as `UMAP`, `isoform_umap`, or another available reduction.

This skill is for feature visualization only. It does not cover QC, merging, annotation, marker statistics, or trajectory inference.

## Use This Skill When

- the user wants to plot genes on a Seurat embedding
- the user wants to plot isoforms on a Seurat embedding
- the user wants to choose the reduction used for plotting
- the user wants to generate multiple feature plots from the same Seurat object without reloading the RDS every time

## Project Source

Execution source of truth:

- `/work/.pi/scripts/draw_features.r`
- `/work/.pi/scripts/draw_features_helpers.r`

Use the helper script as the source of truth for persistent R sessions.
Use the one-off script only for genuine one-shot plotting.

## Current Scope

This skill currently supports:

1. loading one Seurat RDS
2. optionally selecting an assay
3. selecting a reduction
4. drawing one or more gene or isoform features with `SCP::FeatureDimPlot`
5. saving a PDF and summary table
6. supporting repeated plotting from one persistent R session so the same object stays loaded while the user requests additional plots

## Mandatory Rules

1. Run this R-based skill in the `downstream` conda environment.
2. First resolve the input Seurat RDS path.
3. If the user wants to draw multiple plots from the same object, you must keep one R session open so the RDS is loaded only once.
4. In a repeated plotting session, do not restart the script between successive plots unless the user wants to switch to a different RDS.
5. If the user wants a specific assay, confirm it before plotting.
6. If the user wants a specific reduction, confirm it before plotting.
7. Before execution, show the exact final command and wait for confirmation.
8. Do not claim completion until the expected PDF and summary outputs exist and are non-empty.
9. For multiple plots from the same object, repeated shell loops or repeated `Rscript /work/.pi/scripts/draw_features.r ...` calls are not allowed.
10. For multiple plots from the same object, use one persistent `R` session until the user says they are done plotting.
11. If the user first provides one batch of genes or isoforms, draw that batch and then keep the same `R` session open.
12. After finishing one batch, explicitly ask the user whether they want to draw additional genes or isoforms.
13. If the user provides additional genes or isoforms, continue plotting in the same `R` session.
14. Only exit the `R` session after the user explicitly says they do not want to draw anything else.
15. Even if the user initially provides only one gene or one isoform, treat the request as a persistent plotting session unless the user explicitly says they want only a single one-off plot.
16. In a repeated plotting session, one-shot R invocations such as `R --vanilla -f some_script.R` are also not allowed, because they terminate the R process after one scripted run.
17. Reuse the already confirmed workflow-level results root and place default outputs under `14_draw_features` unless the user explicitly asks for a different subdirectory name.

## Entry Procedure

At the start of a feature-drawing request:

1. Resolve the input Seurat RDS path.
2. Decide whether this is:
   - a one-off plot request
   - or a multi-plot session on the same object
3. If it is a multi-plot session, use one persistent R session instead of repeatedly calling `Rscript`.
4. If it is a one-off plot request, confirm:
   - `features`
   - optional `assay`
   - `reduction`
   - `output_file`

If the user initially provides a batch of features and may continue exploring after that first batch, still treat it as a multi-plot session and keep the `R` session alive after the first round of plotting.
If the user initially provides only one gene or one isoform, still treat it as a multi-plot session by default unless the user explicitly says they want a one-off plot only.

## Preferred Mode

### Persistent R Session

Use this mode when the user is likely to ask for several rounds of plots from the same object.

Why:

- the Seurat object is loaded once
- the agent can keep one R terminal alive
- after each plot, the agent can ask the user:
  - which gene or isoform or feature batch to draw next
  - which assay to use
  - which reduction to use
  - whether to continue
- the user can stop the plotting session by saying they do not want more plots, then the agent exits the R session
- repeated `Rscript` calls are not allowed in this situation
- one-shot `R --vanilla -f script.R` execution is also not allowed in this situation

Preferred command pattern:

```bash
conda run -n downstream R
```

In this mode:

- open one R terminal
- load the needed libraries
- `source("/work/.pi/scripts/draw_features_helpers.r")`
- read the Seurat RDS once, preferably with `read_feature_plot_object()`
- wait for the user to provide the next gene or isoform, assay, and reduction
- use `plot_feature_batch()` to draw the requested batch
- use `save_feature_plot()` to save the plot
- ask whether the user wants to draw additional genes or isoforms
- if the user provides another batch, continue in the same R session
- if the user is done, quit the R terminal with `q()`

Do not reload the Seurat object for every plot in this mode.
Do not replace this workflow with a shell `for` loop over repeated `Rscript` invocations.
Do not replace this workflow with a one-shot `R --vanilla -f script.R` execution either.

### One-Off Mode

Use this mode when the user only wants one plot batch.

Preferred command pattern:

```bash
conda run -n downstream Rscript /work/.pi/scripts/draw_features.r \
  --input_rds=/abs/path/to/seurat_object.rds \
  --features=Sox9,Neurog3,Fev,Rbp4 \
  --output_file=/abs/path/to/feature_plot.pdf \
  --assay=RNA \
  --reduction=UMAP
```

## Expected Inputs

Typical inputs for this skill may include:

- `input_rds`
- `features`
- optional `assay`
- `reduction`
- either `output_file` or `output_dir`

## Expected Deliverables

For a normal successful run, produce:

- one PDF per plot request
- one summary TSV per plot request

## Response Style

Be operational and interactive:

- say which Seurat object is being used
- say which assay and reduction are being used
- if the user wants multiple plots, explicitly propose a persistent R session
- in a persistent R session, keep the same session alive and ask whether the user wants to continue after each plot
- if multiple plots are requested from the same object, do not switch to repeated `Rscript` calls
- if the user first gives one batch and then adds more genes or isoforms, keep using the same `R` session until they explicitly say to stop
- if the user starts with only one gene or isoform, still keep the same `R` session alive after drawing it unless they explicitly asked for a one-off plot
- for repeated plotting sessions, do not fall back to `R --vanilla -f ...` or any other one-shot scripted R execution
