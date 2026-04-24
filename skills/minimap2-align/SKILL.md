---
name: minimap2-align
description: Runs project-specific minimap2 splice-aware alignment for long-read RNA-seq, primarily for ONT single-cell reads after BLAZE, with explicit preset selection, tag-preservation handling, and command confirmation.
---

# minimap2 Align

Use this skill for the genome-alignment step in the project workflow.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- long-read RNA reads need splice-aware genome alignment
- downstream processing requires SAM or BAM
- ONT BLAZE-demultiplexed reads need mapping to the genome

## Mandatory Rules

1. Ask minimap2-specific questions only when the workflow reaches this step.
2. Confirm the alignment preset explicitly, typically:
   - `-x splice` for ONT
   - `-x splice:hq` for PacBio CCS if alignment is explicitly requested
3. Confirm thread count explicitly.
4. If BLAZE wrote `CB:Z` and `UB:Z` tags into FASTQ comments, preserve them with `-y`.
5. Derive output SAM under the workflow-level results root.
6. Run the minimap2 command inside a persistent `tmux` session.
7. Before execution, show the exact minimap2 command, the planned `tmux` session name, and wait for confirmation.
8. After launch, do not exit that `tmux` session until the minimap2 command has finished.
9. Execute the confirmed minimap2 command directly in the shell or in `tmux`; do not generate a temporary `.sh` wrapper unless the user explicitly asks for one.
10. Resolve the minimap2 output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
11. Use `02_minimap2` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Default Project Guidance

- ONT standard route: `minimap2 -ax splice -y -t <threads> ref.fa reads.fastq.gz > aln.sam`
- Do not surface low-level scoring and indexing options unless needed.
