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
5. If the user has multiple samples, explicitly remind them to choose a multi-sample mode rather than a one-sample-per-run route.
6. Use fixed conda environments by tool:
   - `BLAZE` -> `blaze`
   - `minimap2` -> `minimap2`
   - `samtools` -> `samtools`
   - `FLAIR` -> `flair`
   - `IsoQuant` -> `isoquant`
7. Run every long-running upstream step inside a persistent `tmux` session.
8. Before executing any step command, show the exact final command, the planned `tmux` session name, and wait for confirmation.
9. After launch, do not exit that `tmux` session until the active step command has finished.
10. Execute confirmed step commands directly in the shell or in `tmux`; do not generate temporary `.sh` wrapper scripts unless the user explicitly asks for them.
11. Confirm one top-level results root once at the start of the upstream workflow, then derive step-specific subdirectories under that root instead of repeatedly asking for new output directories.
12. Show each resolved step output path to the user for confirmation before execution.
13. If dependencies are missing or invalid, stop and report clearly.
14. Do not silently rerun long jobs if outputs may already exist.

## Branches

- `ONT`: `BLAZE -> minimap2 -> samtools -> one of {IsoQuant, bambu, FLAIR}`
- `PacBio`: `one of {IsoQuant, bambu}` from BAM

## Tool Skills To Load

- `../blaze-demultiplex/SKILL.md`
- `../minimap2-align/SKILL.md`
- `../samtools-bam-processing/SKILL.md`
- `../isoquant-quantification/SKILL.md`
- `../bambu-quantification/SKILL.md`
- `../flair-quantification/SKILL.md`

Read the relevant tool skill only when the workflow reaches that step.

## Environment Routing

Use these fixed environment assignments throughout the upstream workflow:

- `blaze-demultiplex` -> `blaze`
- `minimap2-align` -> `minimap2`
- `samtools-bam-processing` -> `samtools`
- `flair-quantification` -> `flair`
- `isoquant-quantification` -> `isoquant`

Do not ask the user to choose an upstream conda environment for these steps unless troubleshooting requires an explicit override.
