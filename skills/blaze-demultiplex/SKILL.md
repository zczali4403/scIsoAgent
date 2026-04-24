---
name: blaze-demultiplex
description: Runs project-specific BLAZE demultiplexing for ONT long-read single-cell RNA-seq, including barcode and UMI detection, whitelist generation, chemistry confirmation, and demultiplexed FASTQ production.
---

# BLAZE Demultiplex

Use this skill for the ONT barcode and UMI preprocessing step.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- the platform is `ONT`
- input is raw long-read FASTQ or FASTQ.gz
- barcode and UMI calling are still needed
- a whitelist and demultiplexed FASTQ are required for downstream steps

## Mandatory Rules

1. Do not use this for PacBio direct-quantification runs.
2. Do not silently rely on BLAZE default chemistry.
3. Confirm one chemistry route explicitly:
   - `--kit-version`
   - or custom barcode route with `--full-bc-whitelist` and `--umi-len`
4. Confirm one whitelist-generation mode explicitly unless running demux-only mode.
5. Derive output paths under the workflow-level results root.
6. Run the BLAZE command inside a persistent `tmux` session.
7. Before execution, show the exact BLAZE command, the planned `tmux` session name, and wait for confirmation.
8. After launch, do not exit that `tmux` session until the BLAZE command has finished.
9. Execute the confirmed BLAZE command directly in the shell or in `tmux`; do not generate a temporary `.sh` wrapper unless the user explicitly asks for one.
10. Resolve the BLAZE output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
11. Use `01_blaze` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Key Decisions To Ask

- chemistry route
- whitelist mode
- thread count
- whether to reuse or overwrite existing outputs

Do not surface low-level tuning unless the user asks or troubleshooting requires it.
