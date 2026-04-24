---
name: isoquant-quantification
description: Runs project-specific IsoQuant single-cell transcript discovery and quantification for ONT aligned-BAM input or PacBio CCS raw FASTQ input, with explicit run organization, barcode route, and command confirmation.
---

# IsoQuant Quantification

Use this skill when the selected Step 4 tool is `IsoQuant`.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- the user selected `IsoQuant`
- ONT input is a sorted and indexed BAM from the upstream alignment route
- or PacBio CCS input is raw FASTQ for direct single-cell processing

## Mandatory Rules

1. Ask IsoQuant-specific questions only when the workflow reaches Step 4.
2. First ask how samples should be organized:
   - one sample per run
   - or multi-sample mode
3. Confirm the correct input route and data type explicitly.
4. Confirm the barcode route explicitly.
5. Confirm thread count explicitly.
6. Derive the output directory under the workflow-level results root.
7. Run every IsoQuant command inside a persistent `tmux` session.
8. Before execution, show the full final command or command block, the planned `tmux` session name, and wait for confirmation.
9. After launch, do not exit that `tmux` session until the active IsoQuant command or command block has finished.
10. Execute confirmed IsoQuant commands directly in the shell or in `tmux`; do not generate temporary `.sh` wrappers unless the user explicitly asks for them.
11. Resolve the IsoQuant output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
12. Use `04_isoquant` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
12. If the user chooses multi-sample mode, do not run one joint single-cell quantification across all samples as the final quantification step.
13. In multi-sample mode, use this required two-stage workflow:
   - first run joint transcript discovery in pseudo-bulk mode on the combined BAM files from all samples with `--mode bulk` to produce one comprehensive `transcript_models.gtf`
   - then run IsoQuant separately for each sample using that GTF as `--genedb` together with `--no_model_construction`
14. In the per-sample quantification stage of multi-sample mode, keep samples separate to avoid barcode collisions while preserving consistent transcript definitions across samples.
15. In the per-sample quantification stage of multi-sample mode, ask the user which single-cell `--mode` should be used for per-sample quantification.

## Project Defaults

- ONT route: `--bam`, `--data_type nanopore`, `--mode tenX_v3`, `--barcoded_bam`
- PacBio route: `--fastq`, `--data_type pacbio_ccs`, `--mode tenX_v3`, plus explicit barcode route

## Multi-Sample Mode Rule

If the user selects multi-sample mode, follow this project-specific pattern:

1. combine BAM inputs across samples for joint transcript discovery in pseudo-bulk mode
   - use `--mode bulk`
2. generate one shared `transcript_models.gtf`
3. rerun IsoQuant separately for each sample using:
   - the sample's own input data
   - `--genedb <transcript_models.gtf>`
   - `--no_model_construction`
4. ask the user which single-cell `--mode` should be used for those per-sample quantification runs

Do not describe the multi-sample route as a single final joint single-cell quantification run.

Do not guess resume, force, or overwrite behavior.
