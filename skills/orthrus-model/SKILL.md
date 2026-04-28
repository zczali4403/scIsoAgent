---
name: orthrus-model
description: Uses Orthrus transcript-sequence embeddings plus pretrained regression and classification heads to predict isoform-level functional properties for one gene at a time in a persistent tmux shell session.
---

# Orthrus Model

Use this skill when the user wants sequence-based isoform functional inference from transcript FASTA sequences using Orthrus embeddings and pretrained task heads.

This skill covers:

1. resolving or preparing a transcript FASTA
2. selecting one gene at a time using the resource-setup transcript-to-gene mapping table
3. extracting a gene-specific transcript FASTA
4. generating Orthrus embeddings for the selected gene's transcripts
5. predicting:
   - mean ribosome load
   - RNA half-life
   - translation efficiency
   - protein localization
   - RNA lifecycle
   - RNA localization
6. aggregating all task outputs into one summary table and one PDF figure
7. continuing gene-by-gene in the same tmux session until the user explicitly says to stop

This skill does not cover upstream quantification itself, Seurat-based downstream analysis, DTU, DeepLocRNA localization prediction, or TranslationAI TIS/TTS prediction.

## Mandatory Session Rule

All work in this skill must be performed in one persistent `tmux` shell session.

Do not split the workflow across detached one-off commands outside `tmux`.
Do not exit the `tmux` session after one gene unless the user explicitly says they are done.

## Environment Routing

Use one fixed conda environment for the whole skill:

- all terminal work, embedding inference, task prediction, and PDF figure generation: `orthrus`

## Project Source Pattern

This skill depends on:

- transcript FASTA preparation in the shell
- one persistent `tmux` shell session
- the local helper script:
  - `../../scripts/orthrus_model_helpers.py`

Prefer using the helper script instead of repeatedly rewriting FASTA extraction, embedding, prediction, and aggregation code inline.

## Use This Skill When

- the user wants Orthrus-based sequence embeddings for isoforms
- the user wants sequence-based functional inference for isoforms from one gene
- the user wants pretrained task-head predictions across multiple transcript properties
- the user wants one PDF figure that consolidates all six Orthrus task outputs

## Entry Procedure

At the start of an Orthrus-model request:

1. Confirm the top-level workflow results root and derive a default output directory under it.
2. Use `15_orthrus_model` as the default subdirectory name unless the user explicitly asks for a different name.
3. Confirm whether the upstream quantification route was:
   - `flair`
   - or `isoquant`
4. If the route was `flair`, ask the user to confirm or provide the transcript FASTA path.
5. If the route was `isoquant`, ask the user to confirm or provide:
   - the transcript annotation GTF path
   - the reference genome FASTA path
   - the output transcript FASTA path
6. Confirm the `resource-setup` dictionary path:
   - `isoform_gene_symbol_dict.csv`
7. Confirm either:
   - one task-manifest CSV/TSV that lists all pretrained Orthrus task heads
   - or the six individual pretrained head paths plus any required label-map JSON files
8. Start or reuse one persistent `tmux` shell session in the `orthrus` environment.
9. Run any FASTA-preparation step first if needed.
10. Ask the user which gene should be analyzed first.
11. Keep the same `tmux` session alive while the user iteratively requests additional genes.

## Workflow

### 1. Resolve Or Prepare Transcript FASTA

If the upstream route was `flair`:

- confirm the existing transcript FASTA path

If the upstream route was `isoquant`:

- prepare a transcript FASTA from the IsoQuant output GTF and the reference genome FASTA

Preferred command pattern:

```bash
conda run -n orthrus gffread annotation.gtf \
  -g genome.fa \
  -w transcripts.fa \
  --w-nocds
```

Rules:

- ask the user to confirm `annotation.gtf`
- ask the user to confirm `genome.fa`
- ask the user to confirm `transcripts.fa`
- do not assume the output path
- show the exact final command before execution
- run this command inside `tmux`

### 2. Extract A Gene-Specific Transcript FASTA

After the full transcript FASTA exists:

- ask the user which gene should be analyzed
- read the `resource-setup` mapping table
- resolve all matching transcript IDs for that gene
- extract those transcripts from the full transcript FASTA
- write a new FASTA containing only that gene's transcripts

Preferred command pattern:

```bash
conda run -n orthrus python ../../scripts/orthrus_model_helpers.py extract-gene-fasta \
  --mapping-csv isoform_gene_symbol_dict.csv \
  --input-fasta transcripts.fa \
  --gene-query GENE_NAME \
  --output-fasta gene_transcripts.fa
```

Rules:

- ask the user to confirm the mapping-table path
- ask the user to confirm the input transcript FASTA path
- ask the user to confirm the gene query
- ask the user to confirm the output gene-specific FASTA path
- show the exact final command before execution
- run this command inside `tmux`

### 3. Infer Orthrus Embeddings

Once the gene-specific FASTA exists, generate one embedding per transcript.

Preferred command pattern:

```bash
conda run -n orthrus python ../../scripts/orthrus_model_helpers.py infer-embeddings \
  --fasta gene_transcripts.fa \
  --output gene_embeddings.npz \
  --csv-output gene_embeddings.csv
```

Rules:

- ask the user to confirm the input gene-specific FASTA path
- ask the user to confirm the output embedding NPZ path
- ask the user whether they also want the optional embedding CSV path
- show the exact final command before execution
- run this command inside `tmux`

### 4. Predict All Six Orthrus Tasks

After embeddings exist:

- run all six pretrained Orthrus task heads on the gene embedding matrix
- write one prediction CSV per task under one confirmed per-gene output directory

Preferred command pattern:

```bash
conda run -n orthrus python ../../scripts/orthrus_model_helpers.py predict-tasks \
  --embedding-npz gene_embeddings.npz \
  --task-manifest orthrus_task_manifest.csv \
  --output-dir task_predictions/
```

Expected tasks include:

- `mrl-sugimoto`
- `rnahl-human`
- `translation-efficiency-human`
- `prot-loc`
- `rna-lifecycle-ietswaart`
- `rna-loc-fazal`

Rules:

- ask the user to confirm the embedding NPZ path
- ask the user to confirm the task-manifest path, or explicitly confirm the six individual probe paths before creating the manifest
- ask the user to confirm the prediction output directory
- show the exact final command before execution
- run this command inside `tmux`

### 5. Aggregate Results And Render A PDF Figure

After all six prediction CSVs exist:

- aggregate all task outputs into one transcript-level summary table
- write a categorical-code map JSON if needed
- render the aggregated plot as a PDF figure while also exporting the aggregated tables

Preferred command pattern:

```bash
conda run -n orthrus python ../../scripts/orthrus_model_helpers.py aggregate-results \
  --result-dir task_predictions/ \
  --output-table aggregated_predictions.csv \
  --output-category-map aggregated_category_code_map.json \
  --output-pdf orthrus_predictions.pdf
```

Rules:

- ask the user to confirm the prediction result directory
- ask the user to confirm the aggregated CSV path
- ask the user to confirm the category-map JSON path
- ask the user to confirm the PDF figure output path
- show the exact final command before execution
- run this command inside `tmux`

### 6. Continue Gene-By-Gene Exploration

After one gene summary is produced:

- ask the user whether they want to analyze another gene
- if yes, keep using the same `tmux` session
- if no, exit the `tmux` session

## Mandatory Rules

1. Run all steps in `tmux`.
2. Use the `orthrus` conda environment for every step in this skill.
3. If the upstream route was `flair`, do not regenerate a transcript FASTA if the user already has a correct one.
4. If the upstream route was `isoquant`, use `gffread` to generate the transcript FASTA before any Orthrus step.
5. Do not assume paths for the transcript FASTA, mapping table, task manifest, output directories, or PDF figure; confirm them.
6. For each selected gene, first extract a gene-specific FASTA, then generate embeddings, then run all six tasks, then aggregate and plot.
7. Before any substantial command is executed, show the exact intended command and wait for confirmation.
8. Do not claim completion until the expected FASTA, embedding, prediction, and PDF artifacts exist and are non-empty.
9. Keep the same `tmux` session alive until the user explicitly says they are done with Orthrus analysis.

## Expected Inputs

Typical inputs for this skill may include:

- upstream quantification route
- transcript FASTA path
- or IsoQuant transcript GTF path plus reference genome FASTA path
- `isoform_gene_symbol_dict.csv`
- one or more user-selected genes
- one task-manifest CSV/TSV or six individual pretrained probe paths
- optional label-map JSON files for classification or multilabel tasks
- output directory

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- transcript FASTA generated from IsoQuant
- one or more gene-specific subset FASTA files
- one or more embedding NPZ files
- optional embedding CSV files
- six per-task prediction CSVs for each analyzed gene
- one aggregated transcript-level prediction table
- one categorical-code map JSON
- one PDF prediction figure per analyzed gene

## Response Style

Be operational and session-aware:

- say whether the workflow is in the FASTA-preparation, embedding, prediction, or aggregation phase
- say which transcript FASTA and mapping table are being used
- say which task-manifest or task-head set is being used
- ask the user which gene they want to inspect next
- keep the same `tmux` shell session alive until the user explicitly says to stop
