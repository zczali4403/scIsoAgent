---
name: summary
description: Summarizes completed long-read single-cell RNA-seq runs using deep results-style integrative writing, with emphasis on post-annotation isoform analysis, sequence-model outputs, transcript-level events, gene-function interpretation, and literature-grounded condition-specific consequences.
---

# Analysis Summary

Use this skill when the user wants a final writeup from completed outputs.

This skill is for reading saved artifacts and writing a manuscript-style, results-focused integrative summary. It is not for rerunning analysis, changing thresholds, or inventing unsupported findings.

## Goal

Default writing target:

- `deep results-style integrative writing`

Default emphasis:

1. annotation and later downstream outputs
2. isoform-focused events
3. sequence-model evidence
4. earlier technical stages only as supporting context

Default formal artifact:

- `16_summary/summary_report.md`

## Use This Skill When

- the user wants a final summary of a completed run
- the user wants manuscript-style results writing with integrated interpretation
- the user wants concrete transcript-level examples explaining isoform shifts
- the user wants sequence-model outputs interpreted together with gene function and condition context

## Core Workflow

1. Confirm the results root and scan readable outputs.
2. Detect completed stages and readable evidence files.
3. Prioritize transcript events supported by post-annotation outputs.
4. Rank prioritized events as:
   - `anchor event`
   - `strong supporting event`
   - `secondary supporting event`
5. For each prioritized event:
   - identify the concrete transcripts involved
   - identify which transcript increases and which decreases across the relevant condition or cell type
   - compare Orthrus, TIS/TTS, and localization outputs across those transcripts
   - classify the cross-model signal as:
     - `convergent`
     - `partially convergent`
     - `discordant`
   - infer the most plausible transcript-state contrast
   - interpret how that contrast may affect the known function of the gene in the relevant condition
6. Perform targeted literature search for prioritized genes, pathways, cell states, or mechanisms.
7. Insert citations directly into the relevant result paragraphs and also provide a final reference list.
8. Draft the final text in event-centered paragraphs, then run adversarial self-review.

## Evidence Rules

- Prefer direct evidence from CSV, TSV, TXT, PDF, PNG, JPG, and readable manifest or log files.
- Do not claim a stage succeeded unless non-empty artifacts exist.
- Do not fabricate findings from filenames or directory names alone.
- Every major biological claim must include:
  - run-specific evidence
  - literature support when available
  - an explicit evidence grade:
    - `A`: direct quantitative support from readable tables
    - `B`: figure-supported trend
    - `C`: cross-step inference supported by at least two independent artifacts or by user-provided prior context plus readable outputs
- If a requested stage has no readable outputs, say so explicitly.

## Sequence-Model Semantics

Use these semantics when interpreting sequence-model outputs.

### Orthrus

#### Mean Ribosome Load (MRL)
- Meaning:
  - predicts the mean number of ribosomes associated with a transcript
- Use:
  - compare relative translational loading tendency across transcripts
- Do not overclaim:
  - do not equate directly with measured protein abundance

#### RNA Half-life
- Meaning:
  - reflects transcript stability-related behavior
- Use:
  - compare whether one transcript appears more stability-compatible than another
- Do not overclaim:
  - do not state a definitive in vivo half-life without direct measurement

#### Translation Efficiency
- Meaning:
  - reflects predicted translation efficiency from ribosome profiling and RNA-seq integration
- Use:
  - compare whether transcripts differ in translation-competent state
- Do not overclaim:
  - do not state that a transcript is experimentally proven to translate more

#### Protein Localization
- Meaning:
  - predicts the subcellular compartments of the encoded protein
- Use:
  - discuss whether transcript switches may alter the functional context of the protein product
- Do not overclaim:
  - do not confuse this with RNA localization

#### RNA Lifecycle
- Meaning:
  - reflects dominant transcript state across release, export, and polysome-loading flow
- Use:
  - discuss whether transcripts appear more associated with export, retention, or translation-associated flow states
- Do not overclaim:
  - treat this as transcript-state tendency, not direct kinetic proof

#### RNA Localization
- Meaning:
  - predicts subcellular compartments where the mRNA is found
- Use:
  - discuss compartmental routing or localization-compatible fate
- Do not overclaim:
  - do not treat as direct localization proof

### DeepLocRNA / transcript localization

- Meaning:
  - predicts transcript-level subcellular localization tendency
- Use:
  - compare whether isoforms from the same gene are predicted to occupy different RNA compartments
- Do not overclaim:
  - do not equate prediction with direct localization measurement

### TranslationAI / TIS-TTS

- Meaning:
  - predicts coding-compatible signals such as ORF presence and initiation/termination-related structure
- Use:
  - discuss whether transcripts differ in coding compatibility or translation-related structure
- Do not overclaim:
  - ORF detection does not prove stable protein production
  - predicted TIS/TTS differences are supportive signals, not proteomic confirmation

## Interpretation Rules

- Prefer within-gene transcript-to-transcript comparison over absolute interpretation.
- Treat sequence-model outputs as interpretable evidence layers, not decorative annotations.
- Do not list Orthrus, TIS/TTS, and localization outputs separately without synthesis.
- For each prioritized event, explicitly answer:
  - what changed
  - why the direction of the change may matter biologically
  - whether the models are convergent, partially convergent, or discordant
  - what transcript state they collectively imply
  - how that state may affect the known function of the gene
  - what that may mean in the compared condition or cell-type context
- If a result is counterintuitive relative to transcript structure, IR status, or DTU direction, do not smooth it away; explain the tension.
- If sequence-model outputs do not converge, say so explicitly.
- Use concrete transcript-level examples rather than only high-level biological themes.
- Do not replace concrete transcript explanation with abstract phrases such as `transcript-state redistribution` unless the transcript example has already been established.
- Prefer direct but still evidence-bounded language.
- Use language such as:
  - `indicates`
  - `points to`
  - `suggests`
  - `is associated with`
  - `shows a more productive / less productive transcript profile`
- Avoid repetitive template phrasing such as:
  - `compatible with`
  - `consistent with`
  - `transcript state` repeated in every paragraph
- Avoid language such as:
  - `proves`
  - `demonstrates`
  - `confirms protein production`

## Gene-Function And Condition Interpretation

For each anchor or strong supporting event:

- explain the known function of the gene
- explain why that function matters in the relevant condition, cell type, or disease context
- explain how the transcript shift could plausibly change that role
- explicitly state what the event may mean for the biological condition being compared
- do not stop at transcript-level interpretation or gene-function description alone
- after describing the gene's known function, explain how the transcript shift could plausibly modify that function in the disease state or cell-type context
- the final sentence of each anchor-event paragraph should state the likely significance of the event for the relevant condition, disease state, or cell-type program
- distinguish explicitly between:
  - disease-state interpretation
  - cell-type-specific interpretation
- do not force a disease interpretation onto a cell-type-specific event, or a cell-identity interpretation onto a within-cell-type disease comparison
- for each anchor or strong supporting event, explicitly answer:
  - what transcript shift occurred
  - what the sequence-model contrasts imply
  - how this could alter the known function of the gene
  - why that may matter for the compared condition or cell-type context
- discuss whether the shift is more indicative of:
  - adaptive remodeling
  - impaired productive output
  - altered RNA handling
  - altered localization-compatible fate
  - altered stress-response behavior

Use literature both for:

- transcript-event interpretation
- gene-function and disease-context interpretation

## Writing Rules

- Write the main body as continuous manuscript-style prose.
- Keep one paragraph for one message.
- The first sentence of each paragraph should state the paragraph's main claim or role.
- Do not organize the final main text as a stepwise recap of module outputs.
- It is acceptable to organize the main text around prioritized transcript events, but each paragraph must contain explicit sequence-model interpretation.
- Insert citations directly into the relevant paragraphs, not only in a terminal bibliography.
- Prefer natural manuscript prose over repeated template expressions.
- Do not give equal narrative weight to all genes.
- Anchor the main narrative on the strongest events.
- Do not use report-style headers such as technical recap or evidence inventory inside the main results body.

## Output Shape

Prefer this order:

1. anchor event paragraph(s)
2. strong supporting event paragraph(s)
3. secondary supporting event paragraph(s)
4. brief integrative synthesis paragraph
5. caveats integrated with claim strength
6. references

Within each event paragraph, prefer:

1. name the event and the direction of transcript change
2. compare key model outputs across the transcripts
3. classify convergence
4. interpret the biological meaning of the directional contrast
5. connect the event to gene function and condition context
6. insert relevant literature citations in the paragraph

## Deep Mode

This skill should use `deep` interpretation mode by default.

Deep mode requires:

- targeted literature search
- event ranking by evidential strength
- transcript-level sequence-model comparison
- gene-function and condition-specific interpretation
- explicit treatment of model convergence or discordance
- substantial manuscript-style prose rather than short recap blocks

## Self-Review

Before finalizing, check:

- is each paragraph centered on one event or one integrative message
- are the strongest events clearly prioritized
- are sequence-model outputs synthesized rather than listed
- is gene function actually used to explain why the event matters
- are any claims stronger than the evidence allows
- are citations inserted where the claim is made
- is any paragraph generic background instead of interpretation of the current run
