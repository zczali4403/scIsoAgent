---
name: samtools-bam-processing
description: Runs project-specific samtools sorting, indexing, quick structural validation, and optional flagstat summary to prepare long-read RNA alignments for downstream quantification.
---

# samtools BAM Processing

Use this skill after alignment output exists and downstream tools require a sorted and indexed BAM.

Source document:

- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/skills/tools/samtools_bam_processing.md`

## Use This Skill When

- minimap2 or another aligner has produced SAM or unsorted BAM
- downstream quantification expects sorted and indexed BAM
- minimal BAM integrity checks or a flagstat summary are needed

## Mandatory Rules

1. Ask samtools-specific questions only when this step is reached.
2. Prefer `samtools sort` directly on aligner output unless the user explicitly wants intermediate BAMs.
3. Confirm thread count for sorting explicitly.
4. Derive sorted BAM, index, and optional QC outputs under the workflow-level results root.
5. Run the samtools command block inside a persistent `tmux` session.
6. If using a coupled block such as `sort` plus `index`, show the full command block and the planned `tmux` session name before execution.
7. Wait for user confirmation before running any samtools command.
8. After launch, do not exit that `tmux` session until the samtools command block has finished.
9. Execute the confirmed samtools command block directly in the shell or in `tmux`; do not generate a temporary `.sh` wrapper unless the user explicitly asks for one.
10. Resolve the samtools output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
11. Use `03_samtools` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Standard Pattern

```bash
samtools sort -@ <threads> -o sample.sorted.bam sample.sam
samtools index sample.sorted.bam
```
