---
name: lr-scrna-upstream-pipeline
description: Orchestrates project-specific upstream long-read single-cell RNA-seq analysis with explicit ONT vs PacBio branching, step ordering, dependency checks, resume rules, and command confirmation before execution.
---

# Long-Read scRNA Upstream Pipeline

Use this skill when the user wants workflow-level orchestration for upstream long-read single-cell RNA-seq analysis in this project.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- the user wants upstream processing from raw long-read reads
- the task requires ONT vs PacBio branch selection
- the task needs step ordering and dependency gating
- the user is not asking for only one isolated tool step

## Mandatory Rules

1. First confirm the platform explicitly: `ONT` or `PacBio`.
2. First confirm one top-level results root name under `results`.
3. Ask questions progressively, only when each step is reached.
4. At Step 4, explicitly ask which compatible tool to use.
5. Run every long-running upstream step inside a persistent `tmux` session.
6. Before executing any step command, show the exact final command, the planned `tmux` session name, and wait for confirmation.
7. After launch, do not exit that `tmux` session until the active step command has finished.
8. Execute confirmed step commands directly in the shell or in `tmux`; do not generate temporary `.sh` wrapper scripts unless the user explicitly asks for them.
9. Confirm one top-level results root once at the start of the upstream workflow, then derive step-specific subdirectories under that root instead of repeatedly asking for new output directories.
10. Show each resolved step output path to the user for confirmation before execution.
11. If dependencies are missing or invalid, stop and report clearly.
12. Do not silently rerun long jobs if outputs may already exist.

## Branches

- `ONT`: `BLAZE -> minimap2 -> samtools -> one of {IsoQuant, bambu, FLAIR}`
- `PacBio`: `one of {IsoQuant, bambu}` from raw FASTQ

## Tool Skills To Load

- `../blaze-demultiplex/SKILL.md`
- `../minimap2-align/SKILL.md`
- `../samtools-bam-processing/SKILL.md`
- `../isoquant-quantification/SKILL.md`
- `../bambu-quantification/SKILL.md`
- `../flair-quantification/SKILL.md`

Read the relevant tool skill only when the workflow reaches that step.
