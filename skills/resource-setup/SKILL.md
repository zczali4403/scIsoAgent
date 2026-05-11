---
name: resource-setup
description: Prepares core downstream resource files for long-read single-cell RNA-seq analysis in this project, especially isoform-to-gene symbol dictionaries and gene-count matrices converted from ENSG IDs to gene symbols.
---

# Resource Setup

Use this skill when the downstream workflow needs standardized resource files before QC, annotation, marker discovery, or isoform-level interpretation.

This skill covers project-specific resource preparation only. It is not the place for downstream QC, clustering, annotation, trajectory analysis, or DTU.

## Use This Skill When

- the user needs downstream resource file creation
- the user needs an isoform-to-gene dictionary
- the user needs gene-level count matrices converted from `ENSGID` to `geneSymbol`
- the downstream workflow needs harmonized feature naming before later modules

## Project Source

Execution source of truth:

- `/work/scIsoAgent/.pi/scripts/resource_file_creation.r`

Prefer the script as the execution source of truth.

## Core Outputs

This skill currently focuses on creating:

1. an isoform-gene-symbol dictionary with:
   - `transcript_id`
   - `gene_id`
   - `gene_symbol`
2. one merged multi-sample gene count matrix rewritten from `ENSGID` rows to `geneSymbol` rows

These outputs are intended to support later downstream modules and to make gene- and isoform-level interpretation easier.

The merged geneSymbol matrix produced here is intended to feed the downstream sample-splitting step before empty-droplet removal.

## Mandatory Rules

1. First confirm all required input paths before drafting commands or code.
2. Prefer calling the project script directly:
   - `/work/scIsoAgent/.pi/scripts/resource_file_creation.r`
3. Collect the script parameters explicitly from the user before execution.
4. First confirm which upstream tool produced the inputs:
   - `isoquant`
   - or `flair`
5. Use distinct output parameters for the two generated files:
   - `dict_output_file`
   - `matrix_output_file`
6. Keep the distinction clear between:
   - tool-specific isoform dictionary inputs
   - per-sample single-cell matrix inputs
   - output resource files derived from them
7. Do not silently assume gene symbols are complete in the reference annotation.
8. If `gene_name` is missing for a record, fall back to `gene_id`.
9. Preserve transcript-level identifiers in the isoform dictionary even when gene symbols are added.
10. Prefer discovering the per-sample single-cell matrix files directly from the chosen run directory under `results` instead of asking the user to type every barcode, feature, and matrix path manually.
11. Only ask the user for explicit per-sample barcode, feature, and matrix paths if automatic discovery from the selected results directory fails or is ambiguous.
12. The agent must infer candidate sample names before execution, present them to the user for confirmation, and then pass them explicitly to the script through `--sample_names`.
13. For matrix merging, read each sample separately, use the confirmed `sample_names` as barcode prefixes, then merge the matrices before gene-symbol conversion.
14. To keep downstream sample splitting safe, normalize sample names by replacing underscores `_` with hyphens `-` before using them as barcode prefixes.
15. If confirmed sample names collapse into duplicates after underscore-to-hyphen normalization, stop and ask the user to provide unique `sample_names` explicitly.
16. When converting the merged gene-level count matrix, deduplicate on `gene_symbol` exactly as defined by the implementation.
17. Run this R-based skill in the `downstream` conda environment.
18. Before executing any substantial command or script, show the exact final command and wait for confirmation.
19. Do not claim completion until the expected output files exist and are non-empty.
20. Reuse the already confirmed workflow-level results root and derive the resource-setup output subdirectory automatically.
21. Use `05_resource_setup` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Entry Procedure

At the start of a resource setup request:

1. Confirm the top-level downstream results directory.
2. Confirm which run directory under `results` should be used for discovery.
3. Confirm the exact script parameters that still need to be passed explicitly:
   - `tool`
   - `reference_gtf_file`
   - `dict_output_file`
   - `matrix_output_file`
   - `sample_names`
4. Confirm the tool-specific dictionary input:
   - if `tool=isoquant`, collect `output_gtf_file`
   - if `tool=flair`, collect `flair_counts_tsv`
5. Discover the per-sample matrix inputs from the selected run directory whenever possible.
6. Infer sample names automatically from the discovered files or sample directories at the agent level, normalize underscores `_` to hyphens `-`, and show the final proposed `sample_names` to the user before execution.
7. Confirm that the discovered per-sample inputs are complete and aligned by sample.
8. Confirm that output paths are writable and belong under the intended project results area.
9. Validate that all declared or discovered input files exist before building commands.

Only after those are clear should you enter the first resource setup step.

## Default Workflow

Use this default order:

1. create the isoform-gene-symbol dictionary using the appropriate upstream tool route
2. discover each sample's matrix inputs from the selected results directory
3. infer candidate sample names at the agent level
4. normalize underscores in sample names to hyphens
5. confirm the final `sample_names` with the user
6. read each sample's matrix inputs separately
7. prefix each sample's barcodes with the confirmed sample name
8. merge all sample matrices by `gene_id`
9. convert the merged matrix from `ENSGID` to `geneSymbol`
10. verify output structure and row naming

## Step 1: Create Isoform-Gene-Symbol Dictionary

Goal:

- build a mapping table linking each transcript or isoform to its `gene_id` and `gene_symbol`

Tool-specific inputs:

- common:
  - `tool`
  - `reference_gtf_file`
- if `tool=isoquant`:
  - `output_gtf_file`
- if `tool=flair`:
  - `flair_counts_tsv`

Expected output:

- a CSV dictionary with:
  - `transcript_id`
  - `gene_id`
  - `gene_symbol`

Script output parameter:

- `dict_output_file`

Project rules:

- if `gene_name` is missing after merge, replace it with `gene_id`
- for `isoquant`, derive `transcript_id` and `gene_id` from the output GTF
- for `flair`, read `flair_quantify.counts.tsv`, take the first column as transcript IDs, and derive `gene_id` by splitting on `_` and taking the last segment

## Step 2: Merge Per-Sample Gene Matrices

Goal:

- construct one merged multi-sample gene matrix before gene-symbol conversion

Preferred discovery behavior:

- do not ask the user to enumerate all per-sample matrix file paths in routine runs
- start from the selected run directory under `results`
- discover the required files automatically, then present the resolved paths before execution

Expected inputs after discovery:

- confirmed `sample_names`
- per-sample barcode files
- per-sample feature files
- per-sample matrix files

Expected behavior:

- for `isoquant`, look for per-sample gene-level matrix triplets such as:
  - `*.gene_grouped_barcode_counts.barcodes.tsv`
  - `*.gene_grouped_barcode_counts.features.tsv`
  - `*.gene_grouped_barcode_counts.matrix.mtx`
- for `flair`, look for per-sample recovery outputs such as:
  - `*.sc_barcodes.tsv`
  - `*.sc_genes.tsv`
  - `*.sc_gene_by_cell.mtx`
- infer candidate sample names from file names or sample directories at the agent level
- normalize underscores in sample names to hyphens before using them as barcode prefixes
- show the final proposed `sample_names` to the user for confirmation before execution
- read each sample's barcode, feature, and matrix files separately
- add `sample_name_` as a prefix to every barcode from that sample
- merge all sample matrices by `gene_id`

Project rule:

- if discovery returns multiple competing candidates for the same sample or missing triplets, stop and report the ambiguity clearly
- only fall back to asking the user for explicit per-sample paths when automatic discovery cannot resolve the inputs cleanly
- if explicit sample-name overrides are provided, they must align with the discovered sample set

## Step 3: Convert The Merged Gene Matrix From ENSGID To GeneSymbol

Goal:

- rewrite the merged gene-level count matrix so downstream steps can use gene symbols instead of ENSG IDs

Expected inputs:

- the merged per-sample gene matrix produced in Step 2
- the dictionary produced in Step 1

Expected output:

- a converted gene count matrix using `gene_symbol` as row names

Script output parameter:

- `matrix_output_file`

Project rule:

- after merging symbol information, remove duplicated rows based on `gene_symbol`

## Script Interface

Preferred execution pattern:

```bash
conda run -n downstream Rscript /work/scIsoAgent/.pi/scripts/resource_file_creation.r \
  --tool=isoquant \
  --output_gtf_file=/abs/path/to/output.gtf \
  --reference_gtf_file=/abs/path/to/reference.gtf \
  --dict_output_file=/abs/path/to/isoform_gene_dict.csv \
  --sample_names=sample1,sample2 \
  --barcodes_paths=/abs/path/to/sample1.barcodes.tsv,/abs/path/to/sample2.barcodes.tsv \
  --features_paths=/abs/path/to/sample1.features.tsv,/abs/path/to/sample2.features.tsv \
  --matrix_paths=/abs/path/to/sample1.matrix.mtx,/abs/path/to/sample2.matrix.mtx \
  --matrix_output_file=/abs/path/to/geneSymbol_gene_count.csv
```

FLAIR route:

```bash
conda run -n downstream Rscript /work/scIsoAgent/.pi/scripts/resource_file_creation.r \
  --tool=flair \
  --flair_counts_tsv=/abs/path/to/flair_quantify.counts.tsv \
  --reference_gtf_file=/abs/path/to/reference.gtf \
  --dict_output_file=/abs/path/to/isoform_gene_dict.csv \
  --sample_names=sample1,sample2 \
  --barcodes_paths=/abs/path/to/sample1.barcodes.tsv,/abs/path/to/sample2.barcodes.tsv \
  --features_paths=/abs/path/to/sample1.features.tsv,/abs/path/to/sample2.features.tsv \
  --matrix_paths=/abs/path/to/sample1.matrix.mtx,/abs/path/to/sample2.matrix.mtx \
  --matrix_output_file=/abs/path/to/geneSymbol_gene_count.csv
```

Skill-level rule:

- the agent should usually discover the matrix paths from the chosen `results` run directory first, infer and confirm `sample_names`, and only then assemble the final command

The current script implements:

- `make_isoquant_isoform_gene_symbol_dict`
- `make_flair_isoform_gene_symbol_dict`
- `make_gene_count_matrix_from_mtx`
- `merge_sample_gene_count_matrices`
- `convert_ENSGID_to_geneSymbol`

Prefer calling the script directly unless the user explicitly asks for in-session code adaptation.

## Expected Inputs

Typical inputs for this skill may include:

- `tool`
- output or extended annotation GTF for IsoQuant
- `flair_quantify.counts.tsv` for FLAIR
- reference annotation GTF
- one selected run directory under `results` used for automatic discovery
- one barcode TSV per sample after discovery
- one feature TSV per sample after discovery
- one matrix MTX per sample after discovery
- optional explicit sample name per sample

## Expected Deliverables

For a normal successful run, produce:

- isoform-gene-symbol dictionary CSV
- merged and converted gene-symbol gene count matrix
- any additional validated resource files explicitly requested by the user

## Response Style

Be operational and explicit:

- say which resource setup step you are in
- say which exact input files are being used
- say which output files are expected next
- show exact commands or code before execution
- do not introduce extra downstream analysis steps inside this skill
