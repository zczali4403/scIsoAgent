---
name: sashimi-plot
description: Draws gene-level sashimi plots from per-sample BAM files after downstream annotation by extracting group-specific cell barcodes from a Seurat object, subsetting BAMs by barcode tag, building ggsashimi input tables, and iteratively plotting one gene at a time in tmux.
---

# Sashimi Plot

Use this skill when the user wants to draw sashimi plots for one gene at a time after downstream annotation, comparing splice-junction usage across cell types, conditions, or condition-restricted cell types.

This skill is for sashimi-plot preparation and plotting only. It does not cover upstream alignment, isoform quantification itself, ordinary feature plots, DTU testing, or isoform classification.

## Mandatory Session Rule

All work in this skill must be performed inside `tmux`.

This skill has two execution phases:

1. Seurat-object inspection and barcode extraction in the `downstream` conda environment
2. BAM subsetting plus `ggsashimi` plotting in the `ggsashimi` conda environment

Do not run substantial commands outside `tmux`.
Do not exit the `tmux` workflow after one gene unless the user explicitly says they are done drawing sashimi plots.

## Environment Routing

Use fixed conda environments by phase:

- Seurat metadata inspection and barcode extraction: `downstream`
- BAM subsetting, manifest creation, and sashimi plotting: `ggsashimi`

If a `tmux` step is entered with `conda activate`, first run:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
```

then activate the needed environment for that phase.

## Project Source Pattern

This skill may use helper scripts written under:

- `/work/.pi/scripts`

If repeated helper logic is needed, place helper scripts there instead of keeping large ad hoc code blocks only in the terminal.

Primary plotting entrypoint:

- `/work/scIsoAgent/.pi/scripts/ggsashimi.py`

Primary local helper script:

- `/work/.pi/scripts/sashimi_plot_helpers.r`

## Use This Skill When

- the user wants to compare splice-junction usage between cell types
- the user wants to compare splice-junction usage between conditions or samples
- the user wants to compare splice-junction usage between conditions within one selected cell type
- the user wants a gene-specific sashimi plot using per-sample BAM files plus a transcriptome GTF
- the user wants to continue plotting multiple genes from the same prepared context

## Entry Procedure

At the start of a sashimi-plot request:

1. Confirm the top-level workflow results root and derive a default output directory under it.
2. Use `14_sashimi_plot` as the default subdirectory name unless the user explicitly asks for a different name.
3. Confirm the annotated Seurat-object RDS path.
4. Start one persistent `tmux` session.
5. In that session, initialize conda and activate `downstream` before reading the Seurat object.
6. In the `downstream` phase, prefer:

```r
source("/work/.pi/scripts/sashimi_plot_helpers.r")
```

7. Read the Seurat object once and inspect available metadata columns.
8. Show the user the relevant available groupings in metadata, especially likely fields such as:
   - `celltype`
   - `sctype_db`
   - `sample`
   - `sample_id`
   - `condition`
9. Ask the user what comparison they want, for example:
   - cell type vs cell type
   - condition vs condition
   - condition vs condition within one selected cell type
10. Resolve which metadata columns define:
   - sample identity
   - group identity for plotting
   - optional cell type restriction
11. Extract per-sample barcode lists for the selected comparison.
12. If sample prefixes are present in Seurat barcodes, strip them before BAM subsetting.
13. Ask the user to confirm each sample's BAM path.
14. Ask the user to confirm which BAM tag stores the barcode:
   - default `CB`
   - or another tag explicitly provided by the user
15. Switch to the `ggsashimi` environment in `tmux` for BAM subsetting and plotting.
16. Build:
   - one or more subset BAM files
   - `input_bams.tsv`
   - `palette.txt`
17. Ask the user to confirm the transcriptome GTF path produced after transcript identification and quantification.
18. Ask the user which gene to plot first.
19. Resolve the gene's genomic region from the GTF and confirm it before plotting.
20. Show the exact final `ggsashimi` command and wait for confirmation.
21. After one gene plot is finished, ask whether the user wants to draw another gene and keep the same `tmux` workflow alive if yes.

## Workflow

### 1. Inspect Seurat Metadata

In the `downstream` environment:

- `source("/work/.pi/scripts/sashimi_plot_helpers.r")`
- read the annotated Seurat-object RDS, preferably with `read_sashimi_seurat_object()`
- list metadata columns with `list_metadata_columns()`
- summarize likely values for grouping columns with `summarize_metadata_column()`
- help the user choose the comparison structure

The skill should explicitly help the user see which `celltype`, `condition`, or `sample` labels are available before asking for the final comparison.

### 2. Extract Per-Sample Barcodes

Based on the selected comparison:

- subset cells from the Seurat object
- split selected cells by sample
- export one barcode list per sample and per plotted group as needed

Expected comparison patterns include:

1. cell type vs cell type
2. condition vs condition using all cells per sample
3. condition vs condition within one selected cell type

Rules:

- barcode extraction must be done per sample
- if Seurat cell names contain sample prefixes, remove the prefixes before matching to BAM tags
- do not assume the barcode tag is always `CB`; confirm it

Prefer using:

- `prepare_sashimi_metadata()`
- `filter_sashimi_metadata()`
- `extract_barcodes_by_sample_and_group()`
- `write_barcode_files()`

### 3. Subset BAM Files By Barcode

For each sample:

- use the sample BAM path confirmed by the user
- subset reads by the chosen barcode tag
- generate a BAM subset for each plotted group when needed
- make sure each BAM used for plotting has a matching `.bai`

Rules:

- default barcode tag is `CB`
- if the user says the BAM uses another tag, use that tag instead
- do not proceed until BAM paths exist and are readable

Preferred command pattern for creating a subset BAM:

```bash
samtools view \
  -@ <threads> \
  -b \
  -D <barcode_tag>:<barcode_file.txt> \
  -o <subset_group.bam> \
  <input_sample.bam>
```

Preferred command pattern for creating the `.bai` index of a subset BAM:

```bash
samtools index \
  -@ <threads> \
  <subset_group.bam>
```

Prefer using:

- `prepare_bam_subset_command()`
- `prepare_bam_index_command()`
- `prepare_subset_bam_manifest()`

### 4. Build `input_bams.tsv`

Construct `input_bams.tsv` with exactly 3 columns:

1. sample name
2. BAM path, with `.bai` index available
3. group name

The third column is used to define group coloring.

Prefer using:

- `write_input_bams_tsv()`

### 5. Build `palette.txt`

Construct `palette.txt` with one color per line:

- one line per plotted group
- use hex color codes such as `#1f77b4`

The number of lines must match the number of plotted groups.

Prefer using:

- `default_group_palette()`
- `write_palette_txt()`

### 6. Resolve Plot Region From GTF

Ask the user which gene to plot.

Then:

- search the confirmed transcriptome GTF for that gene
- resolve the genomic interval needed for `-c`
- show the resolved region to the user before plotting

The skill should use the transcriptome GTF from the quantification workflow, for example from `isoquant` or `flair`, not a generic unrelated annotation by default.

Prefer using:

- `read_transcriptome_gtf()`
- `resolve_gene_region()`

### 7. Run `ggsashimi`

Preferred command pattern:

```bash
python /work/scIsoAgent/.pi/scripts/ggsashimi.py \
  -b input_bams.tsv \
  -c chr10:27040584-27048100 \
  -g annotation.gtf \
  -P palette.txt
```

Supported and recommended options include:

- `-b input_bams.tsv`
- `-c <genomic_region>`
- `-g <annotation.gtf>`
- `-M 10`
- `-C 3`
- `-O 3`
- `--shrink`
- `--alpha 0.25`
- `--base-size=20`
- `--ann-height=4`
- `--height=3`
- `--width=18`
- `-P palette.txt`

Before running, show the exact final command and wait for confirmation.

Prefer using:

- `prepare_ggsashimi_command()`

### 8. Continue Gene-By-Gene Plotting

After one gene is finished:

- ask the user whether they want to plot another gene
- if yes, reuse the same Seurat context, grouping choice, BAM paths, barcode tag, and prepared files when still valid
- if only the gene changes, resolve the new region from the same GTF and rerun plotting
- if the comparison definition changes, rebuild the barcode subsets and `input_bams.tsv`
- if the user is done, exit the workflow

## Expected Inputs

Typical inputs for this skill may include:

- annotated Seurat-object RDS
- sample metadata column
- group metadata column
- optional cell type metadata column
- group labels to compare
- one BAM path per sample
- barcode tag, default `CB`
- transcriptome GTF path
- one or more gene names
- output directory

## Expected Deliverables

For a normal successful run, produce:

- barcode text files for each sample/group subset when needed
- subset BAM files and `.bai` indexes used for plotting
- `input_bams.tsv`
- `palette.txt`
- one or more sashimi plot outputs

## Response Style

Be operational and comparison-aware:

- say which Seurat object is being used
- show which metadata columns are available for grouping
- explain how the requested comparison maps to per-sample barcode extraction
- say which BAM tag is being used for barcode matching
- show the exact `ggsashimi` command before execution
- after one gene, ask whether the user wants to continue with another gene
