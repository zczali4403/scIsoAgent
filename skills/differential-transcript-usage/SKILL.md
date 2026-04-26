---
name: differential-transcript-usage
description: Performs differential transcript usage analysis in one persistent tmux R session using a Seurat object with an isoform assay, including DTU FASTA preparation, pseudobulk aggregation by sample and cell type, GTF subsetting, IsoformSwitchAnalyzeR object initialization, filtering, ORF annotation, DEXSeq-based DTU testing, result export, and switch plotting for user-selected genes.
---

# Differential Transcript Usage

Use this skill when the user wants transcript- or isoform-level differential usage analysis across conditions after a merged Seurat object with an isoform assay is already available.

This skill is for DTU analysis only. It does not replace ordinary differential gene expression, annotation, isoform classification, or isoform diversity profiling.

## Mandatory Session Rule

All R-based analysis in this skill must be performed in one persistent interactive `tmux` R session.

Do not split the DTU workflow across repeated `Rscript` calls.
Do not exit the R session after one substep.
Keep the same `tmux` R terminal alive until the user explicitly says to exit this DTU workflow.

## Project Source Pattern

This skill depends on these packages and project scripts:

- `DRIMSeq`
- `DEXSeq`
- `stageR`
- `IsoformSwitchAnalyzeR`
- `rtracklayer`
- `Seurat`
- `../../scripts/IsoformSwitchAnalyzeR_patch.r`
- `../../scripts/differential_transcript_usage_helpers.r`

At the start of the DTU R session, prefer:

```r
library(DRIMSeq)
library(DEXSeq)
library(stageR)
library(IsoformSwitchAnalyzeR)
library(rtracklayer)
library(Seurat)
source("../../scripts/IsoformSwitchAnalyzeR_patch.r")
source("../../scripts/differential_transcript_usage_helpers.r")
```

## Use This Skill When

- the user wants differential transcript usage analysis between conditions
- the user wants to aggregate isoform expression by `sample_id` and `celltype`
- the user wants to compare transcript usage across biological groups rather than only gene-level abundance
- the user wants IsoformSwitchAnalyzeR plus DEXSeq-based DTU analysis
- the user wants switch plots for selected genes

## Entry Procedure

At the start of a DTU request:

1. Use the `downstream` conda environment for the whole DTU workflow.
2. Start one persistent `tmux` R session and keep it open for the whole workflow.
3. Resolve the merged Seurat-object RDS path that contains the isoform assay.
4. Confirm which quantification tool produced the isoform models:
   - `flair`
   - `isoquant`
   - or another route
5. Confirm the transcript FASTA path:
   - for `flair`, ask the user to confirm the FLAIR combined FASTA path
   - for other tools, ask the user to provide a transcript FASTA that contains one sequence per transcript
6. Confirm the transcript annotation GTF path from the upstream quantification output.
7. Derive the output directory for DTU analysis from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
8. Use `14_dtu` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
9. Confirm which metadata fields encode:
   - sample identity, typically `sample_id`
   - biological condition
   - cell type, typically `celltype`
10. Confirm the exact `resource-setup` output mapping-table path that contains `gene_id` and `gene_symbol`.

## Workflow

### 1. Prepare DTU FASTA

This step is done in the shell before the R workflow starts.

For `flair`, prepare a DTU-ready FASTA by removing the last underscore-delimited suffix from transcript headers.

The source text uses this pattern:

```bash
awk '
/^>/ {
    n = split(substr($0, 2), a, "_")
    if (n > 1) {
        header = a[1]
        for (i = 2; i < n; i++) {
            header = header "_" a[i]
        }
        print ">" header
    } else {
        print
    }
    next
}
{ print }
' input.fa > output.rm_last.fa
```

Rules:

- always show the exact FASTA-preparation command before execution
- ask the user to confirm the input FASTA path and output FASTA path
- preferably generate the shell command with `prepare_flair_dtu_fasta_command()`
- if the tool is not `flair`, do not assume this header rewrite is appropriate; ask the user to provide a ready-to-use transcript FASTA

### 2. Load Seurat Object

Inside the persistent R session:

- load all required packages
- `source("../../scripts/IsoformSwitchAnalyzeR_patch.r")`
- `source("../../scripts/differential_transcript_usage_helpers.r")`
- read the Seurat-object RDS once and keep it in memory for the whole DTU session, preferably with `read_dtu_seurat_object()`

### 3. Build Pseudobulk Isoform Matrix

Aggregate isoform expression across grouped cells using:

- assay: `isoform`
- grouping: `sample_id` and `celltype`

The source text uses:

```r
pseudo.seurat.isoforms <- AggregateExpression(
  seu.obj,
  assays = "isoform",
  return.seurat = FALSE,
  group.by = c("sample_id", "celltype")
)
```

Then:

- convert the aggregated matrix to a data frame
- inspect the resulting grouped column names

Prefer `aggregate_isoform_pseudobulk()` for this step.

### 4. Build Design Matrix For Conditions

Construct the design matrix from grouped pseudobulk column names.

The source text uses the second underscore-delimited field as `condition`.

Rules:

- do not assume this grouping rule is always valid; show the parsed grouping logic to the user and ask for confirmation
- confirm which parsed field corresponds to condition before running DTU
- make `condition` a factor

Prefer `build_dtu_design_matrix()` for this step.

### 5. Subset GTF To The Expressed Transcripts

Read the upstream transcript GTF, then:

- harmonize transcript identifiers in the pseudobulk matrix if needed
- intersect expressed transcript IDs with transcript IDs in the GTF
- subset the GTF to those transcripts
- subset the count matrix to those transcripts
- export the intersected GTF for DTU analysis

The source text uses `rtracklayer::import()` and `export()`.

Always ask the user to confirm the source GTF path and the output DTU GTF path.

Prefer `subset_gtf_for_dtu()` for this step.

### 6. Initialize IsoformSwitchAnalyzeR Object

Ask the user to confirm:

- the DTU GTF path from the previous step
- the transcript FASTA path from Step 1

Then initialize the IsoformSwitchAnalyzeR object with:

- `isoformCountMatrix = subset_cts`
- `isoformRepExpression = subset_cts`
- `designMatrix = samps`
- `isoformExonAnnoation = <DTU GTF>`
- `isoformNtFasta = <DTU FASTA>`

The source text uses `importRdata()` and then `summary()`.

Prefer `initialize_switch_list()` for this step.

### 7. Pre-filter Isoforms

Run `preFilter()` and ask the user to confirm these thresholds before execution, while providing the current defaults:

- `geneExpressionCutoff = 30`
- `isoformExpressionCutoff = 20`
- `IFcutoff = 0.05`
- `removeSingleIsoformGenes = TRUE`
- `acceptedGeneBiotype = NULL`

Do not silently use the defaults without user confirmation.

Prefer `prefilter_switch_list()` for this step.

### 8. Add ORF Information And Gene Symbols

Run `analyzeORF()` on the filtered object.

Then ask the user to confirm the exact resource table path produced by `resource-setup`, and use it to add gene symbols to `isoformFeatures`.

Do not assume one default path from memory.
Do not silently reuse an old resource table unless the user confirms that it is the correct one for the current run.

The source text expects a mapping table with at least:

- `gene_id`
- `gene_symbol`

Prefer `add_orf_and_gene_symbols()` for this step.

### 9. Run DTU Testing

Run DEXSeq-based DTU testing with `isoformSwitchTestDEXSeq()`.

The current source defaults are:

- `reduceToSwitchingGenes = FALSE`
- `reduceFurtherToGenesWithConsequencePotential = FALSE`
- `alpha = 0.01`
- `dIFcutoff = 0.25`
- `onlySigIsoforms = FALSE`

Ask the user to confirm these settings before execution.

After DTU testing:

- run `extractSwitchSummary()`
- run `extractTopSwitches()`
- export the relevant `isoformFeatures` table to CSV

Prefer:

- `run_dtu_test()`
- `export_dtu_results()`

### 10. Draw Switch Plot For A Selected Gene

Ask the user which gene should be plotted.

Also ask the user to confirm:

- `condition1`
- `condition2`
- output PDF path

Then use `switchPlot()` to draw the DTU result for that gene.

Prefer `plot_dtu_switch()` for this step.

Do not assume one fixed gene or one fixed pair of conditions.

## Mandatory Rules

1. Keep the whole DTU workflow in one persistent `tmux` R session until the user says to exit.
2. Before the R workflow starts, first confirm whether a shell FASTA-preparation step is needed.
3. At the beginning of the R session, load the required libraries, `source("../../scripts/IsoformSwitchAnalyzeR_patch.r")`, and `source("../../scripts/differential_transcript_usage_helpers.r")`.
4. Do not reload the Seurat object for every subtask.
5. Do not assume the parsed condition field from grouped pseudobulk column names; show it and ask the user to confirm.
6. Do not assume the GTF path, FASTA path, or resource-table path; ask the user to confirm each of them.
7. Do not silently use preFilter or DTU-testing thresholds; show the defaults and ask the user to confirm them.
8. Prefer helper functions from `../../scripts/differential_transcript_usage_helpers.r` for FASTA command generation, pseudobulk aggregation, design-matrix construction, GTF subsetting, switch-list initialization, filtering, gene-symbol injection, DTU testing, result export, and switch plotting instead of rewriting the same code inline.
9. Do not collapse DTU requests into ordinary differential gene expression.
10. Before any substantial command or code block is executed, show the exact intended code or command and wait for confirmation.
11. Do not claim DTU analysis is complete until the expected result tables and requested switch plots exist and are non-empty.

## Expected Inputs

Typical inputs for this skill may include:

- merged Seurat-object RDS with isoform assay
- transcript FASTA path
- transcript GTF path
- resource-setup mapping table path
- metadata columns for sample, condition, and cell type
- DTU output directory
- optional user-confirmed thresholds
- optional user-selected gene and condition pair for plotting

## Expected Deliverables

Depending on the branch the user chooses, produce only the relevant outputs, for example:

- DTU-ready FASTA
- DTU-specific subset GTF
- IsoformSwitchAnalyzeR object summaries
- DTU result tables
- top-switch summary table
- exported `isoformFeatures` CSV
- switch plot PDF for user-selected genes

## Response Style

Be operational and session-aware:

- state that DTU analysis is being done in one persistent `tmux` R session
- say which Seurat object, GTF, FASTA, and resource table are being used
- say how condition labels are being parsed and ask the user to confirm
- keep the session alive after one gene plot so the user can continue exploring additional genes
