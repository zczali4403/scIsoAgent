---
name: transcript-localization
description: Predicts transcript subcellular localization with DeepLocRNA from transcript FASTA, then summarizes localization classes and gene-specific transcript localization scores in one persistent tmux R session.
---

# Transcript Localization

Use this skill when the user wants transcript-level subcellular localization prediction and follow-up interpretation after upstream long-read transcript quantification.

This skill covers:

1. preparing or resolving a transcript FASTA
2. running DeepLocRNA prediction from the terminal
3. summarizing predicted localization classes
4. selecting one gene at a time using the resource-setup transcript-to-gene mapping table
5. extracting a gene-specific transcript FASTA subset from the original transcript FASTA
6. rerunning DeepLocRNA on that gene-specific subset
7. exploring gene-specific localization scores in a persistent tmux R session

This skill does not cover upstream quantification itself, gene-level QC, clustering, annotation, or DTU.

## Mandatory Session Rule

This skill has two execution phases:

1. terminal preparation and DeepLocRNA prediction in `tmux`
2. one persistent interactive `tmux` R session for downstream plotting and gene-by-gene exploration

Do not split the interactive R analysis across repeated `Rscript` calls.
Do not exit the R session after one gene unless the user explicitly says they are done.

## Environment Routing

Use fixed conda environments by phase:

- non-R terminal work, including `gffread` and DeepLocRNA prediction: `DeepLocRNA`
- R-based follow-up analysis: `downstream`

All terminal work for this skill must run inside `tmux`.

## Project Source Pattern

This skill depends on:

- terminal commands for FASTA preparation and DeepLocRNA prediction
- one persistent `tmux` R session for downstream summaries

Primary local helper script:

- `../../scripts/transcript_localization_helpers.r`

At the start of the R analysis session, prefer:

```r
source("../../scripts/transcript_localization_helpers.r")
```

Then call helper functions from that script instead of repeatedly rewriting the same transcript-localization summary code in the terminal.

## Use This Skill When

- the user wants to predict transcript localization with DeepLocRNA
- the user wants transcript-level localization summaries across all predicted transcripts
- the user wants to inspect one gene's transcript localization pattern
- the user wants heatmaps of localization scores for all transcripts from a selected gene

## Entry Procedure

At the start of a transcript-localization request:

1. Confirm the top-level workflow results root and derive a default output directory under it.
2. Use `15_transcript_localization` as the default subdirectory name unless the user explicitly asks for a different name.
3. Confirm whether the upstream quantification route was:
   - `flair`
   - or `isoquant`
4. If the route was `flair`, ask the user to confirm or provide the transcript FASTA path.
5. If the route was `isoquant`, ask the user to confirm or provide:
   - the transcript annotation GTF path
   - the reference genome FASTA path
   - the output transcript FASTA path
6. Confirm the DeepLocRNA prediction output directory.
7. Confirm the `resource-setup` dictionary path:
   - `isoform_gene_symbol_dict.csv`
8. Run any FASTA-preparation step first if needed.
9. Then run DeepLocRNA prediction.
10. After prediction finishes, start one persistent `tmux` R session in the `downstream` environment.
11. In that session, `source("../../scripts/transcript_localization_helpers.r")`.
12. Keep the same R session alive while the user iteratively requests additional genes.

## Workflow

### 1. Resolve Or Prepare Transcript FASTA

If the upstream route was `flair`:

- confirm the existing transcript FASTA path

If the upstream route was `isoquant`:

- prepare a transcript FASTA from the IsoQuant output GTF and reference genome FASTA

Preferred command pattern:

```bash
conda run -n DeepLocRNA gffread annotation.gtf \
  -g genome.fa \
  -w transcripts.fa \
  --w-nocds
```

Rules:

- ask the user to confirm `annotation.gtf`
- ask the user to confirm `genome.fa`
- ask the user to confirm `transcripts.fa`
- do not assume the output path
- run this command inside `tmux`

### 2. Run DeepLocRNA Prediction

Once the transcript FASTA exists, run DeepLocRNA.

Preferred command pattern:

```bash
conda run -n DeepLocRNA python /opt/DeepLocRNA-main/DeepLocRNA/fine_tuning_deeprbploc_allRNA_prediction.py \
  --fasta transcripts.fa \
  --rna_types mRNA \
  --species Human \
  --output_dir output/
```

Rules:

- ask the user to confirm the transcript FASTA path
- ask the user to confirm the prediction output directory
- show the exact final command before execution
- run the command inside `tmux`

### 3. Start The Persistent R Session

After the DeepLocRNA output exists:

- start one persistent `tmux` R session in the `downstream` environment
- `source("../../scripts/transcript_localization_helpers.r")`
- read the DeepLocRNA prediction table once
- read the `isoform_gene_symbol_dict.csv` mapping table once

Prefer using:

- `read_deeplocrna_predictions()`
- `read_transcript_gene_mapping()`

### 4. Plot Global Prediction-Class Counts

First summarize how many transcripts fall into each `Prediction` category across the full DeepLocRNA output.

Generate:

- one bar plot with:
  - x-axis: prediction class
  - y-axis: transcript count

Prefer using:

- `summarize_prediction_counts()`
- `plot_prediction_class_counts()`

### 5. Gene-Specific Transcript Localization Heatmap

Ask the user which gene to inspect.

Then:

- use the resource-setup mapping table to find all transcripts for that gene
- write the selected transcript IDs to a plain-text file
- extract those transcripts from the original transcript FASTA into a new subset FASTA
- rerun DeepLocRNA on the subset FASTA with `--plot True`
- read the subset DeepLocRNA output for that gene
- build a score matrix across localization compartments from the subset prediction output
- draw a heatmap where:
  - rows are transcripts
  - columns are localization compartments
  - values are DeepLocRNA scores

Preferred extraction-and-rerun pattern:

```bash
conda run -n DeepLocRNA bash -lc 'awk "
BEGIN {
  while ((getline < \"transcript_ids.txt\") > 0) ids[$1] = 1
}
/^>/ {
  header = substr($0, 2)
  keep = (header in ids)
}
keep { print }
" transcripts.fa > transcripts_sub.fa'

conda run -n DeepLocRNA python /opt/DeepLocRNA-main/DeepLocRNA/fine_tuning_deeprbploc_allRNA_prediction.py \
  --fasta transcripts_sub.fa \
  --rna_types mRNA \
  --species Human \
  --output_dir sub_output/ \
  --plot True
```

Prefer using:

- `resolve_gene_transcripts()`
- `write_transcript_ids_file()`
- `prepare_gene_localization_scores()`
- `plot_gene_localization_heatmap()`

### 6. Continue Gene-By-Gene Exploration

After one gene heatmap is produced:

- ask the user whether they want to analyze another gene
- if yes, keep using the same R session
- if no, exit the R session and then exit the `tmux` session

## Mandatory Rules

1. Run all terminal steps in `tmux`.
2. Use the `DeepLocRNA` conda environment for all non-R steps in this skill.
3. Use the `downstream` conda environment for the persistent R phase.
4. If the upstream route was `flair`, do not regenerate a transcript FASTA if the user already has a correct one.
5. If the upstream route was `isoquant`, use `gffread` to generate the transcript FASTA before running DeepLocRNA.
6. Do not assume paths for the transcript FASTA, prediction output, or mapping table; confirm them.
7. At the beginning of the R phase, first `source("../../scripts/transcript_localization_helpers.r")`.
8. Prefer helper functions from `transcript_localization_helpers.r` instead of rewriting the same plotting and filtering code inline.
9. Do not reload the DeepLocRNA output or the mapping table for every gene unless the user changes inputs.
10. For each selected gene, derive a transcript-ID list from the mapping table, extract a subset FASTA from the original transcript FASTA, and rerun DeepLocRNA on that subset before summarizing gene-specific scores.
11. Do not assume one fixed gene; ask the user.
12. Before any substantial command or code block is executed, show the exact intended command or code and wait for confirmation.
13. Do not claim completion until the expected prediction output and requested plots exist and are non-empty.

## Expected Inputs

Typical inputs for this skill may include:

- upstream quantification route
- transcript FASTA path
- or IsoQuant transcript GTF path plus reference genome FASTA path
- DeepLocRNA output directory
- DeepLocRNA prediction table
- `isoform_gene_symbol_dict.csv`
- one or more user-selected genes
- output directory

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- transcript FASTA generated from IsoQuant
- DeepLocRNA prediction output directory
- prediction-class count bar plot
- one or more gene-specific transcript-ID text files
- one or more gene-specific subset FASTA files
- one or more gene-specific DeepLocRNA subset output directories
- gene-specific transcript localization heatmap
- optional per-gene transcript localization score table

## Response Style

Be operational and session-aware:

- state when the workflow is in the terminal DeepLocRNA phase versus the persistent R analysis phase
- say which transcript FASTA and prediction table are being used
- say which mapping table is being used
- ask the user which gene they want to inspect next
- keep the same `tmux` R session alive until the user explicitly says to stop
