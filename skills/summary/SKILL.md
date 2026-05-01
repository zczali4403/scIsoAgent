---
name: summary
description: Summarizes completed long-read single-cell RNA-seq analysis runs with emphasis on post-annotation isoform analysis and sequence-model outputs, then produces an evidence-grounded result summary plus optional discussion-style interpretation based mainly on those downstream findings.
---

# Analysis Summary

Use this skill when the user wants a final summary after one or more upstream or downstream workflow steps have already finished.

This skill is for reading existing outputs and synthesizing them into a concise report. It is not for rerunning analysis, changing thresholds, or inventing results that are not supported by the saved artifacts.

By default, this skill should emphasize results generated after annotation, especially:

- isoform-focused downstream analysis
- sequence-based isoform functional inference
- cross-step interpretation connecting annotation-aware isoform findings to sequence-model outputs

Preferred report route for this skill:

- generate a structured report body from the detected evidence
- produce a Markdown report as the default formal deliverable


## Use This Skill When

- the user wants a summary of a completed workflow run
- the user wants a stepwise recap of what each finished module produced
- the user wants plots, CSV files, TSV files, and other saved outputs interpreted together
- the user wants language that may help write a results section or discussion section
- the user wants a high-level biological interpretation after technical analysis is complete
- the user wants a detailed, literature-grounded biological interpretation
- the user wants deep discussion-style writing supported by online literature search when available

## Scope

This summary skill may read outputs from any completed stage, but the default synthesis priority is:

1. annotation and post-annotation downstream outputs
2. isoform-focused analysis outputs
3. sequence-model outputs
4. earlier technical stages only as supporting context

Readable stages may include:

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
14. transcript localization
15. TIS/TTS analysis
16. Orthrus-based sequence functional prediction

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
15. The default report language for this skill is English unless the user explicitly requests another language.
16. Use `16_summary` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
17. For biological interpretation statements, assign an explicit evidence grade:
   - `A`: direct quantitative evidence from readable tables (CSV/TSV/TXT)
   - `B`: figure-supported trend from readable plots
   - `C`: cross-step inference supported by at least two independent artifacts
17. Every major biological interpretation must include an `Evidence:` line listing concrete source files used.
18. Do not present `B` or `C` evidence as if it were `A`; keep uncertainty explicit.
19. If the user requests `deep` interpretation depth, biological discussion, manuscript-style discussion support, or a detailed summary, perform targeted online literature search before writing the final biological interpretation, unless network tools are unavailable.
20. For `deep` summaries, the report must be substantially more detailed than a stepwise recap and must include explicit biological interpretation plus discussion-style synthesis.
21. When online literature search is available, use it to gather field context, supporting studies, and relevant caveats for the main biological themes emerging from the detected outputs.
22. The final deep discussion must integrate:
   - direct run-specific evidence
   - literature-based background or plausibility support
   - limitations and uncertainty
   - and clearly labeled hypotheses when applicable.
23. When literature search is used, the final report must include an explicit reference list citing the main articles, reviews, or online sources used for biological interpretation.
24. Distinguish clearly between:
   - direct run-specific evidence from the summarized outputs
   - literature-based background or plausibility support
   - and new hypotheses generated by combining the two.
25. If online literature search is unavailable, state that limitation explicitly and restrict the report to artifact-grounded interpretation plus any user-provided literature.
26. When the user requests a detailed summary, prioritize interpretive depth over brevity.
27. If discussion-style interpretation is requested in deep mode, the discussion section should be detailed and substantial; by default, target at least 500 Chinese characters or an equivalently substantial English paragraph block unless the user asks for shorter text.
28. By default, discussion-style writing should be based mainly on annotation-afterward results, especially isoform-focused analysis and sequence-model outputs, rather than on the full workflow.
29. Earlier upstream or QC-related outputs may be mentioned briefly as workflow context, but they should not dominate the discussion unless the user explicitly asks for a full-run discussion.

## Entry Procedure

At the start of a summary request:

1. Confirm the workflow run or top-level results root to summarize.
2. Confirm whether the user wants:
   - technical step summary only
   - biological interpretation only
   - or both technical summary plus biological interpretation
3. Confirm the report language; default to English unless the user explicitly requests another language.
4. Ask whether the user also wants discussion-style writing support.
5. Confirm interpretation depth:
   - `basic`
   - `standard`
   - `deep`
5. If the user requests detailed summary, biological interpretation, or discussion-style writing, default to `deep` interpretation unless the user explicitly asks for a lighter mode.
6. Confirm interpretation boundary:
   - strict evidence-only interpretation
   - or evidence-grounded hypothesis generation
7. Confirm summary emphasis:
   - annotation and later results only
   - annotation plus selected earlier context
   - or full workflow recap
8. If the user does not specify otherwise, default to `annotation and later results only`.
9. If the user wants `deep` interpretation or discussion-style writing, state whether targeted online literature search will be used.
10. Derive the summary output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path.
11. Scan the selected results root and identify which stage directories are present.
12. Within each detected stage directory, identify readable summary artifacts such as CSV, TSV, TXT, PDF, and image files.
13. For `deep` interpretation when network tools are available, identify 2-4 main biological themes from the detected post-annotation outputs, then perform targeted literature search around those themes before writing the final biological interpretation.
14. Report the detected stage list and the main files that will be used as evidence before writing the final summary.

Only after the available evidence has been identified should this skill produce the final synthesis.

## Default Workflow

Use this default order:

1. detect completed stage directories under the selected results root
2. prioritize annotation and post-annotation stage directories
3. collect readable evidence files for prioritized stages
4. summarize prioritized stages separately
5. add only brief supporting context from earlier stages when needed
6. extract cross-stage themes and major findings, emphasizing isoform and sequence-model outputs
7. if requested, draft discussion-style interpretation grounded mainly in the prioritized outputs
8. if `deep` interpretation is requested and network tools are available, perform targeted literature search on the main biological themes
9. assemble the final report in a structure suitable for Markdown delivery

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

### 5. Sequence-Model Summary

If sequence-model modules are present, summarize items such as:

- transcript localization patterns across isoforms
- TIS/TTS differences across isoforms from the same gene
- Orthrus-based predicted differences in transcript properties
- convergence or disagreement between structural isoform findings and sequence-model outputs

Prefer linking these findings back to annotated cell types, condition comparisons, or isoform classes when the evidence supports it.

### 6. Cross-Stage Synthesis

After stepwise summaries, produce a cross-stage synthesis covering:

- the main technical outcome of the workflow
- the main biological signals supported by the post-annotation outputs
- any repeated theme appearing across annotation, isoform classification, diversity, DTU, transcript localization, TIS/TTS, or Orthrus prediction
- any important caveat caused by missing steps, weak outputs, or incomplete summaries

## How To Read Post-Annotation Results

When summarizing results after annotation, interpret them in the following order.

### 1. Annotation As Biological Context

Treat annotation as the biological background layer for later isoform interpretation.

Focus on:

- which cell types or states were identified
- how confident or well-supported the labels appear from saved plots or marker summaries
- which cell populations are most relevant for later isoform findings

Do not stop at restating labels. Use annotation to define the biological context in which later isoform and sequence-model findings should be interpreted.

### 2. Isoform-Focused Results As The Core Signal

Treat post-annotation isoform analysis as the main evidence layer for biological interpretation.

Focus on:

- whether specific genes show multiple isoforms with distinct patterns
- whether certain isoforms are enriched in particular cell types, clusters, or conditions
- whether structural categories, coding status, diversity metrics, or DTU results converge on the same genes or cell populations

When multiple isoform-focused modules point to the same gene, cell type, or condition contrast, treat that convergence as stronger support than any single isolated result.

### 3. Sequence-Model Results As Functional Support

Treat transcript localization, TIS/TTS, and Orthrus predictions as functional interpretation layers that may help explain isoform differences.

Use them to ask:

- do different isoforms from the same gene show different predicted localization patterns?
- do they show differences in predicted start/stop usage or transcript structural endpoints?
- do Orthrus predictions suggest differences in stability, translation-related properties, lifecycle, or localization?

These outputs should usually be interpreted as functional plausibility support rather than definitive proof.

### 4. Integrate The Three Layers

For the strongest summary and discussion, connect:

1. annotated biological context
2. isoform-level structural or usage differences
3. sequence-model-based functional interpretation

Preferred logic pattern:

- `Context`: which annotated cell type, state, or comparison is involved
- `Isoform evidence`: what structural, abundance, diversity, or DTU difference is observed
- `Functional interpretation`: what localization, TIS/TTS, or Orthrus result may help explain that difference

### 5. What To Prioritize In Discussion

By default, prioritize cases where:

- one gene has multiple isoforms with clearly different downstream signatures
- isoform-focused modules and sequence-model modules point in the same direction
- the same biological theme appears across more than one output type

Do not prioritize:

- isolated technical outputs with no biological connection to later findings
- single weak hints unsupported by other post-annotation evidence
- sequence-model predictions that contradict all other evidence unless explicitly framed as uncertainty or a hypothesis

## Discussion-Style Writing

If the user asks for discussion help, produce a separate section that may include:

- the biological significance of the observed annotated celltype or isoform patterns
- possible mechanistic explanations suggested by post-annotation isoform and sequence-model outputs
- links between structural isoform differences, predicted sequence-level properties, and cell-state or condition differences
- caveats, limitations, and follow-up validation ideas

By default, this section should be built mainly from:

- annotation-aware isoform findings
- isoform classification, diversity, and DTU outputs
- transcript localization, TIS/TTS, and Orthrus outputs

Do not default to a full-workflow discussion unless the user explicitly asks for it.

This section must:

- be explicitly labeled as interpretation
- stay consistent with the detected outputs
- avoid overclaiming causality
- note when claims are tentative or hypothesis-generating
- include literature grounding when literature search was performed
- be detailed rather than schematic in `deep` mode

Use a three-layer interpretation format:

1. `Observed` (what is directly seen in artifacts)
2. `Inferred` (evidence-grounded interpretation)
3. `Hypothesis` (testable, explicitly tentative statement)

When interpretation depth is `basic`, keep mostly `Observed` with minimal `Inferred`.
When interpretation depth is `standard`, include `Observed` + `Inferred`.
When interpretation depth is `deep`, include all three layers with explicit evidence grading.

For `deep` discussion requests, the discussion should be a substantial narrative rather than bullet-only text, and should normally exceed 500 Chinese characters (or comparable English length) unless the user explicitly requests brevity.

## Preferred Report Structure

When the user wants a report-oriented output, prefer the following structure:

1. title and run identifier
2. analyzed results root and detected stages
3. evidence inventory
4. annotation and post-annotation technical summary
5. key isoform and sequence-model findings
6. literature-supported biological interpretation
7. discussion-style interpretation focused on annotation-afterward results
8. references
9. caveats and missing outputs
10. candidate follow-up or validation directions

Preferred report outputs at this stage:

- `summary_report.md`
- `evidence_index.tsv` (key interpretation statements mapped to evidence files and grades)
- optionally, when literature search is used: `literature_context.md` or a literature-supported section embedded in `summary_report.md`
- when literature search is used: a `References` section in `summary_report.md`, and optionally `references.bib` or `references.md`

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
- keep the technical recap ordered by workflow stage
- write the report in English by default unless the user explicitly requests another language
- when the user requests a detailed summary, prefer a long-form, evidence-rich narrative instead of a short recap
- if discussion text is requested, present it as a labeled draft rather than as a definitive conclusion
- if the user asks for a report, organize the response in a structure suitable for direct Markdown delivery
