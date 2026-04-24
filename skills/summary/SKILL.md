---
name: summary
description: Summarizes completed long-read single-cell RNA-seq analysis runs by reading step outputs such as plots, CSV/TSV tables, and key Seurat-derived artifacts, then producing a stepwise result summary plus optional discussion-style biological interpretation.
---

# Analysis Summary

Use this skill when the user wants a final summary after one or more upstream or downstream workflow steps have already finished.

This skill is for reading existing outputs and synthesizing them into a concise report. It is not for rerunning analysis, changing thresholds, or inventing results that are not supported by the saved artifacts.

Preferred report route for this skill:

- generate a structured report body from the detected evidence
- produce a Markdown report as the default formal deliverable


## Use This Skill When

- the user wants a summary of a completed workflow run
- the user wants a stepwise recap of what each finished module produced
- the user wants plots, CSV files, TSV files, and other saved outputs interpreted together
- the user wants language that may help write a results section or discussion section
- the user wants a high-level biological interpretation after technical analysis is complete

## Scope

This summary skill may synthesize outputs from any completed stage, including:

1. upstream processing outputs
2. resource setup
3. split gene count by sample
4. empty-droplet removal
5. QC
6. sample merging and optional batch correction
7. isoform object generation and merging
8. add-iso-to-gene integration
9. annotation
10. feature plotting and isoform exploration
11. isoform classification
12. isoform diversity profiling
13. differential transcript usage

## Evidence Sources

Prefer summarizing from saved outputs that the user can inspect directly, such as:

- PDF plots
- PNG or JPG figures if present
- CSV or TSV summary tables
- TXT reports
- RDS-derived summary files if companion text outputs exist
- manifest files
- log files only when they clarify success, failure, warnings, or sample counts

When both figures and tables exist for the same step, use both:

- figures for qualitative interpretation
- tables for explicit counts, rates, categories, and ranked results

## Mandatory Rules

1. Treat this skill as a read-and-summarize workflow, not as a rerun workflow.
2. First confirm which workflow run or results root should be summarized.
3. Prefer scanning the selected run directory automatically instead of asking the user to enumerate every file manually.
4. Before producing conclusions, explicitly report which step directories and output files were actually found.
5. Distinguish clearly between:
   - outputs that exist and were read
   - outputs that were expected but missing
   - and biological interpretation inferred from the available evidence
6. Do not claim that a workflow step succeeded unless output artifacts for that step exist and are non-empty.
7. Do not fabricate biological findings from plot file names or directory names alone.
8. If an artifact cannot be read directly, summarize only what can be supported from adjacent summary files, manifests, or file naming, and say so explicitly.
9. Keep technical summary and biological interpretation separate.
10. Discussion-style text must be framed as interpretation supported by the available outputs, not as a definitive biological claim beyond the evidence.
11. If the user wants manuscript-ready prose, first summarize the evidence and then provide a clearly labeled draft paragraph.
12. If a requested stage has no readable outputs, say so explicitly instead of inventing a placeholder summary.
13. Reuse the confirmed workflow-level results root and summary scope instead of re-asking for already established context.
14. By default, treat the Markdown report as the main output artifact for this skill.

## Entry Procedure

At the start of a summary request:

1. Confirm the workflow run or top-level results root to summarize.
2. Confirm whether the user wants:
   - technical step summary only
   - biological interpretation only
   - or both technical summary plus biological interpretation
3. Ask whether the user also wants discussion-style writing support.
4. Scan the selected results root and identify which stage directories are present.
5. Within each detected stage directory, identify readable summary artifacts such as CSV, TSV, TXT, PDF, and image files.
6. Report the detected stage list and the main files that will be used as evidence before writing the final summary.

Only after the available evidence has been identified should this skill produce the final synthesis.

## Default Workflow

Use this default order:

1. detect completed stage directories under the selected results root
2. collect readable evidence files for each stage
3. summarize each stage separately
4. extract cross-stage themes and major findings
5. if requested, draft discussion-style interpretation grounded in the detected outputs
6. assemble the final report in a structure suitable for Markdown delivery

## Stage Summary Rules

### 1. Upstream Summary

If upstream outputs are present, summarize items such as:

- whether raw reads were demultiplexed
- whether alignment and BAM preparation outputs exist
- which quantification route was used
- whether transcript or gene matrices were produced

Do not over-interpret upstream steps biologically. Keep this section technical unless downstream evidence supports biological interpretation.

### 2. QC And Sample Processing Summary

If downstream QC-related outputs are present, summarize items such as:

- per-sample matrices generated
- empty-droplet removal outputs and retained-cell summaries
- QC thresholding and singlet retention summaries
- sample merge and batch-correction decisions

Prefer explicit counts and summary tables over qualitative statements.

### 3. Cell Type And Feature Summary

If annotation and exploration outputs are present, summarize items such as:

- annotated cell populations or cluster labels
- major marker-supported interpretations
- notable genes or isoforms highlighted in saved plots

When possible, separate:

- cluster-level technical annotation outcomes
- biologically meaningful cell-state patterns

### 4. Isoform-Focused Summary

If isoform modules are present, summarize items such as:

- isoform integration into the gene-level object
- structural category summaries
- coding versus noncoding patterns
- isoform diversity findings
- celltype-specific isoform findings
- DTU comparisons and switch-like patterns

Prefer quoting categories, counts, rankings, and group comparisons from saved summary tables where available.

### 5. Cross-Stage Synthesis

After stepwise summaries, produce a cross-stage synthesis covering:

- the main technical outcome of the workflow
- the main biological signals supported by the outputs
- any repeated theme appearing across annotation, isoform classification, diversity, or DTU
- any important caveat caused by missing steps, weak outputs, or incomplete summaries

## Discussion-Style Writing

If the user asks for discussion help, produce a separate section that may include:

- the biological significance of the observed celltype or isoform patterns
- possible mechanistic explanations suggested by the outputs
- links between structural isoform differences and cell-state or condition differences
- caveats, limitations, and follow-up validation ideas

This section must:

- be explicitly labeled as interpretation
- stay consistent with the detected outputs
- avoid overclaiming causality
- note when claims are tentative or hypothesis-generating

## Preferred Report Structure

When the user wants a report-oriented output, prefer the following structure:

1. title and run identifier
2. analyzed results root and detected stages
3. evidence inventory
4. stepwise technical summary
5. cross-stage biological interpretation
6. discussion-style interpretation
7. caveats and missing outputs

Preferred report outputs at this stage:

- `summary_report.md`

## Expected Deliverables

For a normal successful summary request, deliver only what the user asked for, for example:

- a detected-stage inventory
- a stepwise technical summary
- a biological interpretation summary
- a discussion-style draft paragraph or short section
- a clear list of missing or unreadable stages that limited interpretation
- `summary_report.md` as the default formal report artifact

## Response Style

Be explicit and evidence-aware:

- say which run directory was summarized
- say which stages were detected
- say which files or file types supported the summary
- separate direct evidence from interpretation
- keep the technical recap concise and ordered by workflow stage
- if discussion text is requested, present it as a labeled draft rather than as a definitive conclusion
- if the user asks for a report, organize the response in a structure suitable for direct Markdown delivery
