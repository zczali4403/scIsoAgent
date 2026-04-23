# IsoAgent User Guide

`IsoAgent` is a long-read single-cell RNA-seq workflow assistant for this project.

It supports:

- upstream processing from raw ONT or PacBio reads
- downstream analysis from quantified outputs and Seurat objects

This repository is the portable workflow bundle used by the project.

## What This Bundle Contains 📦

The repository mainly includes:

- `skills/`: workflow modules
- `scripts/`: analysis scripts and helper scripts

Typical layout:

```text
README.md
scripts/
skills/
prompt_example.txt
```

Keep this structure unchanged.

## What IsoAgent Can Do

IsoAgent currently supports:

- ONT long-read single-cell RNA-seq
- PacBio CCS long-read single-cell RNA-seq

### Upstream ⬆️

Upstream processing turns raw reads into transcript- or gene-level outputs.

#### ONT

```text
Raw FASTQ
  -> BLAZE
  -> minimap2
  -> samtools
  -> one of:
       - IsoQuant
       - bambu
       - FLAIR
```

#### PacBio

```text
Raw FASTQ
  -> one of:
       - IsoQuant
       - bambu
```

### Downstream ⬇️

Downstream analysis currently includes:

- resource setup
- gene-count splitting by sample
- empty-droplet removal
- QC
- sample merging and optional batch correction
- isoform-object generation and merging
- adding isoform assay into gene-level Seurat objects
- annotation
- feature plotting
- isoform exploration
- isoform classification
- isoform diversity analysis
- differential transcript usage

## How To Use It 🚀

In normal use:

- tell the agent what you want to do
- say whether you are starting from raw data or from existing outputs
- provide the main input paths
- confirm the proposed command before execution

Examples:

- "Run the upstream ONT workflow from these FASTQ files"
- "Start downstream analysis from this merged Seurat object"
- "Run only IsoQuant on these aligned BAM files"
- "Plot these genes on the merged object"

## How Outputs Are Organized 📁

IsoAgent uses one top-level `results` directory for a workflow run.

You confirm that root once, and each step then creates its own subdirectory under it.

Typical default subdirectory names:

- `01_blaze`
- `02_minimap2`
- `03_samtools`
- `04_isoquant`
- `04_flair`
- `04_bambu`
- `05_resource_setup`
- `06_split_gene_count`
- `07_remove_droplets`
- `08_qc_pipeline`
- `09_merge_remove_batch`
- `10_generate_iso_obj`
- `11_merge_iso_obj`
- `12_add_iso_to_gene`
- `13_annotation`
- `14_draw_features`
- `14_explore_isoforms`
- `14_isoform_classification`
- `14_isoform_diversity`
- `14_dtu`

## What You Should Prepare 📝

Depending on the task, you may need some or all of the following:

- raw FASTQ files
- aligned BAM files
- reference genome FASTA
- reference annotation GTF or GFF
- quantified transcript or gene matrices
- Seurat RDS objects
- a conda environment with the required tools installed

For downstream work, it is also helpful to know:

- which results directory you want to use
- which Seurat object should be the main input
- which metadata columns represent sample, condition, cluster, or cell type

## Typical Starting Points 🧭

### Start From Raw Data

If you are starting from raw reads, provide:

- platform: `ONT` or `PacBio`
- input file paths
- reference FASTA
- annotation GTF or GFF
- top-level results directory

### Start From Existing Outputs

If you already completed part of the workflow, provide:

- which outputs already exist
- where they are located
- which downstream step you want next

Common starting points include:

- merged gene count matrix
- per-sample count matrices
- Seurat RDS after droplet removal
- QC-completed Seurat objects
- merged Seurat object with isoform assay

## Portability

If you want to move this workflow bundle to another machine or another project location:

- copy the whole bundle
- keep `skills/` and `scripts/` as sibling directories

Do not copy only `skills/` without `scripts/`.

## Summary ✅

Use this workflow bundle if you want a structured way to run long-read scRNA-seq analysis in this project.

The basic pattern is:

- start from your real input files or existing outputs
- confirm one top-level results directory
- let each workflow step create its own subdirectory
- review and confirm commands before they run

If you are unsure where to begin, describe your current inputs and your goal.
