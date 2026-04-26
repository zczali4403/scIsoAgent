---
name: samtools-bam-processing
description: Runs project-specific samtools sorting, indexing, quick structural validation, and optional flagstat summary to prepare long-read RNA alignments for downstream quantification.
---

# samtools BAM Processing

Use this skill after alignment output exists and downstream tools require a sorted and indexed BAM.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- minimap2 or another aligner has produced SAM or unsorted BAM
- downstream quantification expects sorted and indexed BAM
- minimal BAM integrity checks or a flagstat summary are needed

## Mandatory Rules

1. Use the `samtools` conda environment for this skill's commands.
2. Ask samtools-specific questions only when this step is reached.
3. Prefer `samtools sort` directly on aligner output unless the user explicitly wants intermediate BAMs.
4. Confirm thread count for sorting explicitly.
5. Derive sorted BAM, index, and optional QC outputs under the workflow-level results root.
6. Run the samtools command block inside a persistent `tmux` session.
7. If using a coupled block such as `sort` plus `index`, show the full command block and the planned `tmux` session name before execution.
8. Wait for user confirmation before running any samtools command.
9. After launch, do not exit that `tmux` session until the samtools command block has finished.
10. Execute the confirmed samtools command block directly in the shell or in `tmux`; do not generate a temporary `.sh` wrapper unless the user explicitly asks for one.
11. Resolve the samtools output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
12. Use `03_samtools` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Standard Pattern

```bash
samtools sort -@ <threads> -o sample.sorted.bam sample.sam
samtools index sample.sorted.bam
```
