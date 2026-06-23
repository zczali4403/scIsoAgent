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
6. Reuse the existing aligned and indexed BAM from the workflow rather than defaulting to realignment.
7. Confirm splice-junction support source explicitly:
   - short-read STAR `SJ.out.tab`
   - or long-read-derived junctions from `intronProspector`
8. If the workflow uses `flair combine`, require the command to include `-c` / `--convert_gtf`.
9. In this project, do not treat `-c` as optional once `flair combine` is selected, because the combined transcriptome should also be converted to GTF.
10. Run every FLAIR command block inside a persistent `tmux` session.
11. If the session is entered with `conda activate`, first run `source "$(conda info --base)/etc/profile.d/conda.sh"` inside `tmux`, then `conda activate flair`.
12. Confirm all command blocks before execution, including `intronProspector`, `flair transcriptome`, optional `flair combine`, `flair quantify`, and matrix recovery, and show the planned `tmux` session name.
13. After launch, do not exit that `tmux` session until the active FLAIR command block has finished.
14. Execute confirmed FLAIR commands directly in the shell or in `tmux`; do not generate temporary `.sh` wrappers unless the user explicitly asks for them.
15. Resolve the FLAIR output subdirectory automatically from the already confirmed workflow-level results root and ask the user to confirm the resolved path rather than asking for a new output directory.
16. Use `04_flair` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.
17. After `flair quantify`, run the project script:
   - `/work/scIsoAgent/.pi/scripts/recover_sc_matrix.py`
18. Before matrix recovery, verify read names follow:
   - `barcode_umi#read_uuid_strand`
19. In multi-sample runs, run `/work/scIsoAgent/.pi/scripts/recover_sc_matrix.py` only on the per-sample `isoform.read.map.txt` files produced by `flair quantify`.
20. Do not use the multi-sample `combine` result as the input for `recover_sc_matrix.py`.
21. Because the newer `flair transcriptome`, `flair combine`, and `flair quantify` options may require FLAIR from GitHub master, check the installed FLAIR help if an option is rejected instead of silently switching to the older `correct` / `collapse` route.

## Project Defaults

- Default output directory: `<results_root>/04_flair`
- Default alignment route: reuse existing upstream aligned BAM
- Optional realignment route: direct `minimap2` spliced alignment only when confirmed by the user
- Default long-read junction support: `intronProspector` followed by `--junction_bed <junctions.bed> --junction_support 2`
- Default quantification mode: `flair quantify --stringent --check_splice`
- Multi-sample comparative route:
  1. generate one transcriptome per sample
  2. combine sample transcriptomes with `flair combine -c`
  3. quantify all samples against the combined transcriptome
  4. recover single-cell matrices per sample from per-sample read-map outputs

## Optional Alignment Template

Default project behavior is to reuse the aligned BAM from the upstream route. Only use this template if the user confirms realignment is needed.

Recommended FLAIR-protocol minimap2 command:

```bash
minimap2 -ax splice -s 40 -G 350k -t <threads> --MD --secondary=no \
  <genome_fasta> <sample_reads_fasta_or_fastq> \
| samtools view -hb - \
| samtools sort -@ <threads> -o <sample>.genomealigned.bam -
samtools index <sample>.genomealigned.bam
```

For stranded libraries, add `-uf` to the minimap2 options.

## Splice-Junction Support Templates

### Short-Read STAR Junctions

If the user has a STAR `SJ.out.tab` file, pass it directly to `flair transcriptome`:

```text
--junction_tab <SJ.out.tab>
```

### Long-Read Junctions

Generate supported junctions from the aligned BAM:

```bash
intronProspector \
  --genome-fasta=<genome_fasta> \
  --intron-bed6=<sample>.IPjunctions.bed \
  -C 0.0 \
  <sample_genomealigned_bam>
```

Then pass the junction BED to `flair transcriptome`:

```text
--junction_bed <sample_or_combined.IPjunctions.bed> --junction_support 2
```

For multiple samples, prefer one shared junction file for the dataset when possible. Concatenating per-sample long-read junction BED files is acceptable when a shared file is needed.

## Transcriptome Generation Templates

Use `flair transcriptome` for each sample-specific aligned BAM. This is the preferred new FLAIR route and replaces the older `correct` plus `collapse` workflow.

With short-read junction support:

```bash
flair transcriptome \
  -g <genome_fasta> \
  -f <annotation_gtf> \
  --junction_tab <SJ.out.tab> \
  -b <sample_genomealigned_bam> \
  -o <sample_transcriptome_prefix>
```

With long-read junction support:

```bash
flair transcriptome \
  -g <genome_fasta> \
  -f <annotation_gtf> \
  --junction_bed <IPjunctions.bed> \
  --junction_support 2 \
  -b <sample_genomealigned_bam> \
  -o <sample_transcriptome_prefix>
```

Expected transcriptome files for later steps are typically:

```text
<sample_transcriptome_prefix>.isoforms.fa
<sample_transcriptome_prefix>.isoforms.bed
```

## Multi-Sample Combine Rule

If the user selects multi-sample comparative mode, first create a combine manifest pointing to per-sample transcriptome BED files.

Manifest shape:

```text
<sample1>	<sample1_transcriptome_prefix>.isoforms.bed
<sample2>	<sample2_transcriptome_prefix>.isoforms.bed
<...additional_samples>	<...>.isoforms.bed
```

More stringent combine template:

```bash
flair combine \
  -m <combine_manifest.txt> \
  -o <combined_prefix> \
  -c
```

Less stringent combine template:

```bash
flair combine \
  -m <combine_manifest.txt> \
  -o <combined_prefix> \
  -p 0 \
  -f 1 \
  -s \
  -c
```

Do not show a `flair combine` command for confirmation unless `-c` or `--convert_gtf` is present.

Expected combined transcriptome files are typically:

```text
<combined_prefix>.isoforms.fa
<combined_prefix>.isoforms.bed
<combined_prefix>.isoforms.gtf
```

## Quantification Templates

`flair quantify` uses raw FASTA/FASTQ reads, which may be gzipped. Confirm the raw read path for each sample before quantification.

Read manifest shape:

```text
<sample1>	<sample1_reads.fasta_or_fastq_or_gz>
<sample2>	<sample2_reads.fasta_or_fastq_or_gz>
<...additional_samples>	<...reads.fasta_or_fastq_or_gz>
```

Multi-sample quantification against the combined transcriptome:

```bash
flair quantify \
  -r <reads_manifest.txt> \
  -i <combined_prefix>.isoforms.fa \
  --isoform_bed <combined_prefix>.isoforms.bed \
  --stringent \
  --check_splice \
  -o <quant_prefix>
```

One-sample quantification without `flair combine`:

```bash
flair quantify \
  -r <reads_manifest.txt> \
  -i <sample_transcriptome_prefix>.isoforms.fa \
  --isoform_bed <sample_transcriptome_prefix>.isoforms.bed \
  --stringent \
  --check_splice \
  -o <sample_quant_prefix>
```

Expected per-sample read-map outputs are typically:

```text
<sample>.isoform.read.map.txt
```

Verify the actual output names before matrix recovery.

## Single-Cell Matrix Recovery Template

Before recovery, verify that read names in each per-sample `isoform.read.map.txt` follow:

```text
barcode_umi#read_uuid_strand
```

Run recovery once per sample read-map file:

```bash
python /work/scIsoAgent/.pi/scripts/recover_sc_matrix.py \
  --map-file <sample>.isoform.read.map.txt \
  --out-dir <sample_matrix_output_dir> \
  --prefix <sample>
```

Expected outputs per sample:

```text
<sample>.sc_transcript_by_cell.mtx
<sample>.sc_barcodes.tsv
<sample>.sc_transcripts.tsv
<sample>.sc_gene_by_cell.mtx
<sample>.sc_genes.tsv
```

Do not run matrix recovery on the combined transcriptome output.

## Version Check

If a FLAIR command fails because an option is unavailable, check the installed interface:

```bash
flair --help
flair transcriptome --help
flair combine --help
flair quantify --help
```

If the installed FLAIR version lacks the newer protocol options, stop and report the version mismatch to the user.
