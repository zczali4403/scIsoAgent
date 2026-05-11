---
name: flair-quantification
description: Runs project-specific ONT-only FLAIR transcriptome generation, optional multi-sample combine, quantification, and post-quantification single-cell matrix recovery with recover_sc_matrix.py.
---

# FLAIR Quantification

Use this skill when the selected Step 4 tool is `FLAIR`.

Primary skill file:

- `./SKILL.md`

## Use This Skill When

- the user selected `FLAIR`
- the platform is `ONT`
- aligned BAM and original raw reads are both available
- the goal is isoform discovery plus FLAIR-style quantification

## Mandatory Rules

1. Use the `flair` conda environment for all FLAIR commands in this skill.
2. Do not use this skill for PacBio.
3. Ask FLAIR-specific questions only when the workflow reaches the FLAIR step.
4. First ask how samples should be organized:
   - one sample per run
   - or multi-sample comparative run
5. If the user provides multiple samples or multiple sample-specific input files, explicitly remind them that they should choose the multi-sample comparative run.
6. Reuse the existing aligned BAM from the workflow rather than defaulting to realignment.
7. Confirm splice-junction support source explicitly.
8. If the workflow uses `flair combine`, require the command to include `-c` / `--convert_gtf`.
9. In this project, do not treat `-c` as optional once `flair combine` is selected, because the combined transcriptome should also be converted to GTF.
10. Run every FLAIR command block inside a persistent `tmux` session.
11. If the session is entered with `conda activate`, first run `source "$(conda info --base)/etc/profile.d/conda.sh"` inside `tmux`, then `conda activate flair`.
12. Confirm all command blocks before execution, including `intronProspector`, `flair transcriptome`, optional `flair combine`, `flair quantify`, and matrix recovery, and show the planned `tmux` session name.
13. After launch, do not exit that `tmux` session until the active FLAIR command block has finished.
14. Execute confirmed FLAIR commands directly in the shell or in `tmux`; do not generate temporary `.sh` wrappers unless the user explicitly asks for them.
15. Resolve the FLAIR output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
16. Use `04_flair` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
17. When showing a `flair combine` command for confirmation, explicitly check that `-c` or `--convert_gtf` is present.
18. After `flair quantify`, run the project script:
   - `/work/scIsoAgent/.pi/scripts/recover_sc_matrix.py`
19. Before matrix recovery, verify read names follow:
   - `barcode_umi#read_uuid_strand`
20. In multi-sample runs, run `/work/scIsoAgent/.pi/scripts/recover_sc_matrix.py` only on the per-sample `isoform.read.map.txt` files produced by `flair quantify`.
21. Do not use the multi-sample `combine` result as the input for `recover_sc_matrix.py`.

## Output Handling

Derive transcriptome, quantify, combine, and matrix-recovery outputs under the workflow-level results root.

If `flair combine` is used, expect both the combined transcriptome output and the converted GTF output because `-c` / `--convert_gtf` is required in this project.

For multi-sample runs, matrix recovery is a per-sample post-quantification step tied to each sample's `isoform.read.map.txt` output from `flair quantify`.
