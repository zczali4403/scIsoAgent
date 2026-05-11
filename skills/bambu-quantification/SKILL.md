---
name: bambu-quantification
description: Runs project-specific bambu single-cell transcript discovery and quantification for ONT aligned-BAM input or PacBio raw FASTQ input, with explicit run organization, chemistry handling, and command confirmation.
---

# bambu Quantification

Use this skill when the selected Step 4 tool is `bambu`.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- the user selected `bambu`
- transcript discovery and transcript or gene count matrices are needed
- ONT input is a demultiplexed and aligned BAM
- or PacBio input is raw FASTQ in direct single-cell mode

## Mandatory Rules

1. Ask bambu-specific questions only when the workflow reaches the bambu step.
2. First ask how samples should be organized:
   - one sample per run
   - or multi-sample joint run
3. If the user provides multiple samples or multiple sample-specific input files, explicitly remind them that they should choose the multi-sample joint run.
4. Confirm the correct input route explicitly:
   - `--bams` for ONT aligned-BAM route
   - `--reads` for direct raw-read route
5. For PacBio direct runs, confirm chemistry explicitly.
6. Confirm thread count explicitly.
7. Derive `--outdir` under the workflow-level results root.
8. Run the final Nextflow command or command block inside a persistent `tmux` session.
9. Before execution, show the exact final Nextflow command or command block, the planned `tmux` session name, and wait for confirmation.
10. After launch, do not exit that `tmux` session until the bambu workflow command has finished.
11. Execute the confirmed bambu workflow command directly in the shell or in `tmux`; do not generate a temporary `.sh` wrapper unless the user explicitly asks for one.
12. Resolve the bambu output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
13. Use `04_bambu` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Runtime Notes

- This project assumes bambu runs through Nextflow.
- Container runtime details are operational, not biological decisions, and should not be surfaced unless needed.
