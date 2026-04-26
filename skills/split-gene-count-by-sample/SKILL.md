---
name: split-gene-count-by-sample
description: Splits the merged geneSymbol gene-count matrix produced after resource setup into per-sample matrices by parsing sample names from barcode prefixes before empty-droplet removal.
---

# Split Gene Count By Sample

Use this skill when the downstream workflow needs to split the merged gene-level count matrix into one matrix per sample before empty-droplet removal.

This skill is for per-sample matrix splitting only. It does not create the merged matrix, run empty-droplet removal, or run downstream QC.

## Use This Skill When

- the user has completed `resource-setup`
- the user has one merged geneSymbol gene-count matrix
- barcode names encode sample identity as the substring before the first `_`
- the next step requires one gene-count matrix per sample

## Project Source

Execution source of truth:

- `../../scripts/split_gene_count_by_sample.r`

## Core Outputs

This skill currently produces:

1. one per-sample geneSymbol gene-count matrix for each sample detected in the merged matrix
2. one manifest file listing:
   - sample name
   - barcode count
   - output file path

## Mandatory Rules

1. Prefer calling the project script directly:
   - `../../scripts/split_gene_count_by_sample.r`
2. Prefer detecting the merged geneSymbol matrix from the most recent `resource-setup` output when available.
3. If automatic detection is missing, ambiguous, or clearly wrong, ask the user to provide the merged matrix path manually.
4. Infer sample names from barcode prefixes before the first `_`.
5. If any barcode does not match the expected `sample_barcode` pattern, stop and report the offending columns clearly.
6. Run this R-based skill in the `downstream` conda environment.
7. Before execution, show the exact final command and wait for confirmation.
8. Do not claim completion until the expected per-sample matrices and manifest file exist and are non-empty.
9. Reuse the already confirmed workflow-level results root and derive the split output subdirectory automatically; do not ask the user to type a fresh output directory unless they explicitly want to override it.
10. Use `06_split_gene_count` as the default subdirectory name under the confirmed workflow-level results root unless the user explicitly asks for a different name.

## Entry Procedure

At the start of a split request:

1. Resolve the merged geneSymbol gene-count matrix path:
   - either detect it from the most recent `resource-setup` output
   - or ask the user to provide `gene_count_file`
2. Derive the split output directory from the already confirmed workflow-level results root and ask the user to confirm the resolved path:
   - `output_dir`
3. Ask whether barcode prefixes should be stripped after splitting:
   - recommended default: keep prefixes (`strip_prefix=false`)
4. Validate that the input matrix exists before proposing commands.

## Script Interface

Preferred execution pattern:

```bash
conda run -n downstream Rscript ../../scripts/split_gene_count_by_sample.r \
  --gene_count_file=/abs/path/to/geneSymbol_gene_count.csv \
  --output_dir=/abs/path/to/output_dir \
  --strip_prefix=false
```

Required parameters:

- `gene_count_file`
- `output_dir`

Optional parameters:

- `strip_prefix`

## Expected Deliverables

For a normal successful run, produce:

- one `*_geneSymbol_gene_count.csv` file per sample
- `split_manifest.tsv`

## Response Style

Be operational and explicit:

- say which merged matrix is being split
- say how sample names are being inferred
- say which output directory is being used
- show the exact command before execution
