---
name: lr-scrna-downstream-pipeline
description: Orchestrates project-specific downstream long-read single-cell RNA-seq analysis after upstream quantification, including resource setup, droplet cleanup, QC, sample merging, isoform matrix integration, cell type identification, isoform interpretation, isoform classification, isoform diversity profiling, and differential transcript usage.
---

# Long-Read scRNA Downstream Pipeline

Use this skill when the user wants downstream analysis for long-read single-cell RNA-seq in this project.

This skill is the downstream entry point. It defines scope, stage ordering, handoff rules, and the current module map. It should stay stable even if specific downstream methods, packages, or module-level skills change later.

## Use This Skill When

- the user wants downstream analysis after upstream long-read single-cell processing
- the user already has transcript-by-cell and/or gene-by-cell matrices
- the user wants cell clustering, annotation, isoform interpretation, isoform classification, isoform diversity analysis, or DTU
- the user asks workflow-level questions about downstream analysis order, dependencies, or deliverables

## Current Scope

This downstream pipeline currently includes the following planned modules:

1. resource file creation
2. split gene count matrix by sample
3. quality control
4. sample merging and optional batch correction
5. isoform expression matrix integration
6. cell type identification
7. isoforms of interest exploration
8. isoform classification
9. isoform diversity profiling
10. differential transcript usage
11. final run summary and discussion-style interpretation

This is a living workflow definition. The module list, order, and boundaries may be revised as the project matures.

## Mandatory Rules

1. Treat this skill as a routing and orchestration layer, not as the place for package-specific implementation detail.
2. At the start of a downstream request, first confirm the top-level downstream results directory.
3. Before planning later modules, first confirm which upstream outputs are available and usable.
4. Do not assume all downstream modules are required for every analysis.
5. Ask questions progressively. Do not front-load parameters for modules that the user may not run.
6. If a module-specific skill exists, read that skill before proposing detailed commands or code.
7. If a requested downstream module does not yet have a dedicated skill, say so explicitly and proceed with a minimal project-consistent plan instead of inventing hidden structure.
8. Keep gene-level and isoform-level analyses clearly separated when the distinction matters for interpretation.
9. Before executing any substantial command or script, show the exact final command and wait for confirmation.
10. Do not claim a module is complete until its expected output artifacts exist and are non-empty.

## Entry Procedure

At the start of a new downstream analysis request:

1. Confirm one top-level downstream results root under `results`.
2. After that root is confirmed, derive per-skill and per-step output subdirectories automatically under it instead of repeatedly asking the user to provide new output directories.
3. Show each resolved output path to the user for confirmation before execution.
2. Confirm the upstream source and available inputs, such as:
   - gene-by-cell matrix
   - transcript-by-cell or isoform-by-cell matrix
   - cell metadata
   - transcript or isoform annotation
   - extended annotation if it already exists
   - dimensional reduction, cluster, or preliminary label files if they already exist
3. Confirm the analysis level the user wants:
   - gene-focused
   - isoform-focused
   - or joint gene + isoform downstream analysis
4. Confirm which downstream modules are in scope for this run.
5. Validate that all declared input files exist before drafting module-level commands.

Only after those are clear should the workflow enter the first required downstream module.

## Default Stage Order

Use this default order unless the user explicitly requests a different entry point or already has completed outputs for some stages:

1. resource file creation
2. split gene count matrix by sample
3. quality control
4. sample merging and optional batch correction
5. isoform expression matrix integration
6. cell type identification
7. isoforms of interest exploration
8. isoform classification
9. isoform diversity profiling
10. differential transcript usage
11. final run summary and discussion-style interpretation

This ordering is a default, not a rigid law. Some analyses may start from an intermediate stage if validated prior outputs already exist.

## Module Map

### 1. Resource File Creation

Use this module to prepare the core downstream resource set that later steps depend on.

Typical inputs may include:

- upstream expression matrices
- metadata tables
- transcript-to-gene mapping
- isoform annotation tables
- reference annotation or extended annotation outputs

Typical outputs may include:

- standardized downstream input manifest
- harmonized feature metadata
- harmonized cell metadata
- transcript-to-gene or isoform-to-gene mapping resources

### 2. Split Gene Count Matrix By Sample

Use this module after `resource-setup` and before empty-droplet removal.

Typical goals:

- take the merged geneSymbol gene-count matrix
- infer sample identity from barcode prefixes before the first `_`
- write one per-sample matrix for downstream empty-droplet removal
- establish the per-sample input set that `remove-droplets` should iterate over

Current dispatch rule for this stage:

- if the task is splitting the merged gene matrix into per-sample matrices, read:
  - `../split-gene-count-by-sample/SKILL.md`

### 3. Quality Control

Use this module for early downstream cleanup and quality control.

Current planned QC scope includes:

- empty droplet removal
- routine QC metrics and filtering

Project rule:

- after `split gene count matrix by sample`, run empty-droplet removal separately for each sample matrix rather than once on the merged matrix

Current dispatch rule for this stage:

- if the task is workflow-level downstream QC, read:
  - `../qc-pipeline/SKILL.md`
- if the task is empty droplet removal, read:
  - `../remove-droplets/SKILL.md`

Do not silently assume the exact QC method or threshold set. These may differ across datasets and may later move into additional dedicated QC skills.

### 4. Sample Merging And Optional Batch Correction

Use this module after per-sample QC and doublet removal, and before integrated multi-sample downstream interpretation.

Typical goals:

- merge all per-sample QC-completed singlet Seurat objects
- inspect sample structure in a merged UMAP before integration
- decide whether batch correction is necessary
- if needed, run Harmony-based integration and then save the final merged object

Current dispatch rule for this stage:

- if the task is merging QC-completed per-sample Seurat objects or deciding whether to remove batch effects, read:
  - `../merge-remove-batch/SKILL.md`

### 5. Isoform Expression Matrix Integration

Use this module when isoform-level downstream analysis is needed.

Typical goals:

- add isoform expression matrices into the downstream analysis object
- keep gene-level and isoform-level assays aligned to the same cells
- preserve feature identity and transcript-to-gene relationships
- create one combined object for joint downstream gene-plus-isoform analysis

Current dispatch rule for this stage:

- if the task is creating per-sample isoform Seurat objects from transcript barcode-feature-matrix files and matching them to gene-level QC-retained cells, read:
  - `../generate-iso-obj/SKILL.md`
- if the task is merging per-sample isoform Seurat objects after `generate-iso-obj`, read:
  - `../merge-iso-obj/SKILL.md`
- if the task is adding the merged isoform assay and isoform UMAP back into the merged gene object, read:
  - `../add-iso-to-gene/SKILL.md`

### 6. Identifying Cell Types

Use this module to assign biological meaning to clusters or cell states.

Typical evidence sources may include:

- marker genes
- marker isoforms
- known reference annotations
- project-specific biological context

Current dispatch rule for this stage:

- if the task is cluster-to-celltype annotation on the merged Seurat object, read:
  - `../annotation/SKILL.md`

### 7. Exploring Isoforms of Interest

Use this module when the user wants focused follow-up on selected genes or transcripts.

Typical outputs may include:

- per-cell or per-group isoform expression summaries
- usage patterns across clusters or conditions
- prioritized candidate isoforms for validation

Current dispatch rule for this stage:

- if the task is plotting genes or isoforms on a Seurat reduction such as `UMAP` or `isoform_umap`, read:
  - `../draw-features/SKILL.md`
- if the task is focused exploratory analysis of selected genes or isoforms, including isoform-count summaries, top isoform-rich genes, selected-gene isoform inspection, or IsoVis export preparation, read:
  - `../exploring-isoforms-of-interest/SKILL.md`

### 8. Isoform Classification

Use this module to classify isoforms into biologically or structurally meaningful categories.

Current focus in this project includes:

- SQANTI3-based structural classification
- structural category summaries across cell types
- structural subcategory summaries across cell types
- coding-versus-noncoding proportion summaries across cell types
- celltype-specific isoform detection

Current dispatch rule for this stage:

- if the task is SQANTI-based isoform structural classification plus downstream summary of structural categories, subcategories, coding status, or celltype-specific isoforms, read:
  - `../isoform-classification/SKILL.md`

### 9. Profiling Isoform Diversity

Use this module to quantify and compare isoform complexity across cells, clusters, or conditions.

Examples of possible goals:

- diversity within a gene across cell states
- diversity across clusters
- changes in isoform repertoire complexity

Current dispatch rule for this stage:

- if the task is entropy-based isoform diversity analysis in a selected cell type, including entropy heatmaps, high- and low-entropy gene summaries, or per-cell isoform proportion pie charts, read:
  - `../profiling-isoform-diversity/SKILL.md`

### 10. Differential Transcript Usage

Use this module for transcript- or isoform-level comparative analysis across conditions, groups, or trajectories.

Keep DTU clearly distinct from ordinary differential gene expression. Do not collapse transcript-usage questions into gene-level DE by default.

Current dispatch rule for this stage:

- if the task is IsoformSwitchAnalyzeR- or DEXSeq-based differential transcript usage analysis, including DTU FASTA preparation, pseudobulk aggregation, GTF subsetting, DTU testing, or switch plotting, read:
  - `../differential-transcript-usage/SKILL.md`

### 11. Final Summary

Use this module when the user wants a final synthesis after analysis outputs have been generated.

Typical goals:

- summarize what each completed stage produced
- integrate evidence from plots and summary tables
- highlight the major technical and biological findings
- draft discussion-style interpretation if requested

Current dispatch rule for this stage:

- if the task is summarizing a completed workflow run, interpreting saved outputs across steps, or drafting discussion-style takeaways from existing result files, read:
  - `../summary/SKILL.md`

## Module Dispatch Rules

When a downstream request enters a concrete module:

- first check whether a dedicated downstream module skill already exists
- if it exists, read it before asking detailed questions or drafting commands
- if it does not exist, stay explicit that the project currently has pipeline-level guidance only for that module

Current downstream module skills already available:

- resource file creation:
  - `../resource-setup/SKILL.md`
- split gene count matrix by sample:
  - `../split-gene-count-by-sample/SKILL.md`
- downstream QC workflow:
  - `../qc-pipeline/SKILL.md`
- empty droplet removal:
  - `../remove-droplets/SKILL.md`
- sample merging and optional batch correction:
  - `../merge-remove-batch/SKILL.md`
- isoform object generation:
  - `../generate-iso-obj/SKILL.md`
- isoform object merging:
  - `../merge-iso-obj/SKILL.md`
- add merged isoform assay back to merged gene object:
  - `../add-iso-to-gene/SKILL.md`
- feature plotting on Seurat reductions:
  - `../draw-features/SKILL.md`
- cell type annotation:
  - `../annotation/SKILL.md`
- focused isoform exploration:
  - `../exploring-isoforms-of-interest/SKILL.md`
- isoform classification:
  - `../isoform-classification/SKILL.md`
- isoform diversity profiling:
  - `../profiling-isoform-diversity/SKILL.md`
- differential transcript usage:
  - `../differential-transcript-usage/SKILL.md`
- final run summary and discussion-style interpretation:
  - `../summary/SKILL.md`

## Expected Inputs

Depending on the branch and quantification route, downstream analysis may consume:

- gene-by-cell matrix
- transcript-by-cell or isoform-by-cell matrix
- cell barcode list
- sample metadata
- cluster labels
- low-dimensional embeddings
- transcript annotation GTF/GFF
- extended annotation outputs from upstream tools

Do not assume every run includes all of these.

## Expected Deliverables

For a normal successful downstream run, deliver only the artifacts relevant to the modules actually requested, for example:

- cleaned and standardized downstream resource files
- QC summaries and filtered objects
- merged multi-sample Seurat objects with or without batch correction
- integrated gene-level and isoform-level analysis objects
- cell type labels
- isoform-focused summary tables and plots
- isoform classification outputs
- isoform diversity summaries
- DTU result tables and supporting summaries
- a final cross-step summary and optional discussion-style interpretation

## Response Style

Be operational and explicit:

- say which downstream module you are in
- say which validated inputs are being used
- say what output artifacts are expected next
- show exact commands before execution
- do not over-specify methods that the project has not fixed yet
- keep the pipeline document stable even when module implementation details evolve
