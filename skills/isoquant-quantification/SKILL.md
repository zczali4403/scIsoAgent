---
name: isoquant-quantification
description: Runs project-specific IsoQuant single-cell transcript discovery and quantification for ONT aligned-BAM input or PacBio CCS BAM input, with explicit run organization, barcode route, and command confirmation.
---

# IsoQuant Quantification

Use this skill when the selected Step 4 tool is `IsoQuant`.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- the user selected `IsoQuant`
- ONT input is a sorted and indexed BAM from the upstream alignment route
- or PacBio CCS input is a BAM file

## Mandatory Rules

1. Use the `isoquant` conda environment for all IsoQuant commands in this skill.
2. Ask IsoQuant-specific questions only when the workflow reaches Step 4.
3. First ask how samples should be organized:
   - one sample per run
   - or multi-sample mode
4. If the user provides multiple samples or multiple sample-specific input files, explicitly remind them that they should choose `multi-sample mode`.
5. Confirm the correct input route and data type explicitly.
6. Confirm the barcode route explicitly.
7. For PacBio BAM input, first determine whether barcode and UMI information are already present in BAM tags.
8. If a PacBio BAM does not already contain barcode and UMI tags such as `CB` and `UB`, but the reads still retain the single-cell barcode / UMI / polyT structure, use IsoQuant single-cell barcode calling such as `--mode tenX_v3` together with `--barcode_whitelist`.
9. If a PacBio BAM already contains barcode and UMI tags, do not repeat barcode calling; instead use `--barcoded_bam` so IsoQuant reads barcode and UMI directly from BAM tags.
10. If a barcoded PacBio BAM uses non-default tags instead of `CB` and `UB`, explicitly confirm and pass `--barcode_tag` and `--umi_tag`.
11. Confirm thread count explicitly.
12. Derive the output directory under the workflow-level results root.
13. Run every IsoQuant command inside a persistent `tmux` session.
14. Before execution, show the full final command or command block, the planned `tmux` session name, and wait for confirmation.
15. After launch, do not exit that `tmux` session until the active IsoQuant command or command block has finished.
16. Execute confirmed IsoQuant commands directly in the shell or in `tmux`; do not generate temporary `.sh` wrappers unless the user explicitly asks for them.
17. Resolve the IsoQuant output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
18. Use `04_isoquant` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
19. If the user chooses multi-sample mode, do not run one joint single-cell quantification across all samples as the final quantification step.
20. In multi-sample mode, use this required two-stage workflow:
   - first run joint transcript discovery in pseudo-bulk mode on the combined BAM files from all samples with `--mode bulk` to produce one comprehensive `transcript_models.gtf`
   - then run IsoQuant separately for each sample using that GTF as `--genedb` together with `--no_model_construction`
21. In the per-sample quantification stage of multi-sample mode, keep samples separate to avoid barcode collisions while preserving consistent transcript definitions across samples.
22. In the per-sample quantification stage of multi-sample mode, ask the user which single-cell `--mode` should be used for per-sample quantification.

## Project Defaults

- ONT route: `--bam`, `--data_type nanopore`, `--mode tenX_v3`, `--barcoded_bam`
- PacBio route without prewritten barcode tags: `--bam`, `--data_type pacbio_ccs`, `--mode tenX_v3`, `--barcode_whitelist <whitelist>`
- PacBio route with prewritten barcode tags: `--bam`, `--data_type pacbio_ccs`, `--barcoded_bam`
- PacBio route with non-default barcode tags: add `--barcode_tag <tag>` and `--umi_tag <tag>`

## PacBio Barcode Route

For PacBio data, the BAM structure determines which IsoQuant single-cell route should be used:

1. If the input BAM is an ordinary PacBio HiFi / CCS BAM without single-cell barcode and UMI tags, but the reads still retain the 10x barcode / UMI / polyT structure, use IsoQuant barcode calling directly in single-cell mode, for example `--mode tenX_v3`, and provide `--barcode_whitelist`.
2. If the input BAM already contains barcode and UMI tags written during upstream preprocessing, do not repeat barcode calling. Use `--barcoded_bam` so IsoQuant reads barcode and UMI from BAM tags.
3. IsoQuant defaults to `CB` and `UB`. If the BAM uses other tags such as `CR` / `UR` or `XC` / `XM`, pass them explicitly with `--barcode_tag` and `--umi_tag`.

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
