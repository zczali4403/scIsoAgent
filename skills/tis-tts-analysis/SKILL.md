---
name: tis-tts-analysis
description: Predicts transcript translation initiation and termination sites with TranslationAI from transcript FASTA, then summarizes gene-specific TIS/TTS coordinates in one persistent tmux R session.
---

# TIS TTS Analysis

Use this skill when the user wants transcript-level translation start and stop site prediction with TranslationAI after upstream long-read transcript quantification.

This skill covers:

1. preparing or resolving a transcript FASTA
2. running TranslationAI from the terminal
3. selecting one gene at a time using the resource-setup transcript-to-gene mapping table
4. extracting and displaying the predicted TIS/TTS coordinates for the selected gene's transcripts in one plot

This skill does not cover upstream quantification itself, gene-level QC, clustering, annotation, or DTU.

## Mandatory Session Rule

This skill has two execution phases:

1. terminal preparation and TranslationAI prediction in `tmux`
2. one persistent interactive `tmux` R session for downstream gene-by-gene visualization

Do not split the interactive R analysis across repeated `Rscript` calls.
Do not exit the R session after one gene unless the user explicitly says they are done.

## Environment Routing

Use fixed conda environments by phase:

- non-R terminal work, including `gffread` and TranslationAI prediction: `TranslationAI`
- R-based follow-up analysis: `downstream`

All terminal work for this skill must run inside `tmux`.

## Project Source Pattern

This skill depends on:

- terminal commands for FASTA preparation and TranslationAI prediction
- one persistent `tmux` R session for downstream summaries

Primary local helper script:

- `/work/scIsoAgent/.pi/scripts/tis_tts_analysis_helpers.r`

At the start of the R analysis session, prefer:

```r
source("/work/scIsoAgent/.pi/scripts/tis_tts_analysis_helpers.r")
```

Then call helper functions from that script instead of repeatedly rewriting the same TIS/TTS parsing and plotting code in the terminal.

## Use This Skill When

- the user wants to predict translation initiation and termination sites with TranslationAI
- the user wants to inspect one gene's predicted start and stop positions across its transcripts
- the user wants a single plot showing multiple transcripts from one gene together

## Entry Procedure

At the start of a TIS/TTS analysis request:

1. Confirm the top-level workflow results root and derive a default output directory under it.
2. Use `15_tis_tts_analysis` as the default subdirectory name unless the user explicitly asks for a different name.
3. Confirm whether the upstream quantification route was:
   - `flair`
   - or `isoquant`
4. If the route was `flair`, ask the user to confirm or provide the transcript FASTA path.
5. If the route was `isoquant`, ask the user to confirm or provide:
   - the transcript annotation GTF path
   - the reference genome FASTA path
   - the output transcript FASTA path
6. Confirm the TranslationAI prediction output directory.
7. Confirm the `resource-setup` dictionary path:
   - `isoform_gene_symbol_dict.csv`
8. Run any FASTA-preparation step first if needed.
9. Then run TranslationAI prediction.
10. After prediction finishes, start one persistent `tmux` R session in the `downstream` environment.
11. If the session is entered with `conda activate`, first run `source "$(conda info --base)/etc/profile.d/conda.sh"` inside `tmux`, then `conda activate downstream`, then start `R`.
12. In that session, `source("/work/scIsoAgent/.pi/scripts/tis_tts_analysis_helpers.r")`.
13. Keep the same R session alive while the user iteratively requests additional genes.

## Workflow

### 1. Resolve Or Prepare Transcript FASTA

If the upstream route was `flair`:

- confirm the existing transcript FASTA path

If the upstream route was `isoquant`:

- prepare a transcript FASTA from the IsoQuant output GTF and reference genome FASTA

Preferred command pattern:

```bash
conda run -n TranslationAI gffread annotation.gtf \
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

### 2. Run TranslationAI Prediction

Once the transcript FASTA exists, run TranslationAI.

Preferred command pattern:

```bash
conda run -n TranslationAI translationai \
  -I transcripts.fa \
  -t 0.5,0.5 \
  --plain_fasta \
  --output_dir outputdir
```

Rules:

- ask the user to confirm the transcript FASTA path
- ask the user to confirm the prediction output directory
- show the exact final command before execution
- run the command inside `tmux`

### 3. Start The Persistent R Session

After the TranslationAI ORF output exists:

- start one persistent `tmux` R session in the `downstream` environment
- `source("/work/scIsoAgent/.pi/scripts/tis_tts_analysis_helpers.r")`
- read the TranslationAI ORF output once
- read the `isoform_gene_symbol_dict.csv` mapping table once

Prefer using:

- `read_translationai_orf_output()`
- `read_transcript_gene_mapping()`

### 4. Gene-Specific TIS/TTS Visualization

Ask the user which gene to inspect.

Then:

- use the resource-setup mapping table to find all transcripts for that gene
- subset the TranslationAI output to those transcripts
- extract the predicted TIS and TTS coordinates
- draw one plot containing all matched transcripts from that gene

The plot should show:

- one horizontal track per transcript
- the predicted translation initiation site
- the predicted translation termination site

Prefer using:

- `resolve_gene_transcripts()`
- `prepare_gene_tis_tts_table()`
- `plot_gene_tis_tts()`

### 5. Continue Gene-By-Gene Exploration

After one gene plot is produced:

- ask the user whether they want to analyze another gene
- if yes, keep using the same R session
- if no, exit the R session and then exit the `tmux` session

## Mandatory Rules

1. Run all terminal steps in `tmux`.
2. Use the `TranslationAI` conda environment for all non-R steps in this skill.
3. Use the `downstream` conda environment for the persistent R phase.
4. If the upstream route was `flair`, do not regenerate a transcript FASTA if the user already has a correct one.
5. If the upstream route was `isoquant`, use `gffread` to generate the transcript FASTA before running TranslationAI.
6. Do not assume paths for the transcript FASTA, prediction output, or mapping table; confirm them.
7. At the beginning of the R phase, first `source("/work/scIsoAgent/.pi/scripts/tis_tts_analysis_helpers.r")`.
8. Prefer helper functions from `tis_tts_analysis_helpers.r` instead of rewriting the same parsing and plotting code inline.
9. Do not reload the TranslationAI output or the mapping table for every gene unless the user changes inputs.
10. Do not assume one fixed gene; ask the user.
11. Before any substantial command or code block is executed, show the exact intended command or code and wait for confirmation.
12. Do not claim completion until the expected prediction output and requested plots exist and are non-empty.

## Expected Inputs

Typical inputs for this skill may include:

- upstream quantification route
- transcript FASTA path
- or IsoQuant transcript GTF path plus reference genome FASTA path
- TranslationAI output directory
- TranslationAI ORF output file
- `isoform_gene_symbol_dict.csv`
- one or more user-selected genes
- output directory

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- transcript FASTA generated from IsoQuant
- TranslationAI prediction output directory
- one or more gene-specific TIS/TTS coordinate tables
- one or more gene-specific TIS/TTS plots

## Response Style

Be operational and session-aware:

- state when the workflow is in the terminal TranslationAI phase versus the persistent R analysis phase
- say which transcript FASTA and ORF output file are being used
- say which mapping table is being used
- ask the user which gene they want to inspect next
- keep the same `tmux` R session alive until the user explicitly says to stop
