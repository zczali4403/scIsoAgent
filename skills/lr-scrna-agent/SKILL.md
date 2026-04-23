---
name: lr-scrna-agent
description: Orchestrates the full long-read single-cell RNA-seq workflow for this project, covering both upstream processing and downstream analysis. Use when the user wants a single entry skill for ONT or PacBio long-read single-cell RNA-seq from raw data through downstream interpretation.
---

# Long-Read scRNA Agent

Use this skill as the top-level entry point for the whole long-read single-cell RNA-seq workflow in `/home4/zcz/zcz_proj1/IsoAgent/agent_proj`.

This skill covers both:

- upstream processing
- downstream analysis

It is a routing and orchestration skill. It should decide whether the user is asking for upstream work, downstream work, or a cross-stage continuation, then hand off to the correct specialized entry skill.

## Use This Skill When

- the user wants a complete long-read single-cell RNA-seq workflow
- the user has not clearly separated upstream and downstream steps
- the user wants to continue from upstream outputs into downstream analysis
- the user asks workflow-level questions about the whole project
- the user wants one unified entry point instead of choosing separate upstream or downstream skills

## Project Scope

Supported platforms:

- ONT long-read single-cell RNA-seq
- PacBio CCS long-read single-cell RNA-seq

Supported workflow layers:

1. upstream processing
2. downstream analysis

## Primary Entry Skills

This top-level skill should route to:

- upstream entry:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-upstream-agent/SKILL.md`
- downstream entry:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-downstream-pipeline/SKILL.md`

## Upstream Coverage

Upstream processing includes:

- raw long-read single-cell input handling
- ONT barcode recovery with BLAZE when applicable
- alignment with minimap2
- BAM processing with samtools
- transcript or isoform discovery and quantification with:
  - IsoQuant
  - bambu
  - FLAIR
- transcript-by-cell and gene-by-cell matrix generation

Current supported upstream routes:

- ONT: `BLAZE -> minimap2 -> samtools -> one of {IsoQuant, bambu, FLAIR}`
- PacBio: `one of {IsoQuant, bambu}` from raw FASTQ

## Downstream Coverage

Downstream analysis currently includes:

1. resource file creation
2. split gene count matrix by sample
3. empty-droplet removal and downstream QC
4. sample merging and optional batch correction
5. isoform expression matrix integration
6. cell type identification
7. isoforms of interest exploration
8. isoform classification
9. isoform diversity profiling
10. differential transcript usage

## Mandatory Operating Rules

Follow these rules exactly:

1. First determine whether the user is asking for upstream work, downstream work, or continuation from upstream into downstream.
2. If the request is upstream, immediately route conceptually to:
   - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-upstream-agent/SKILL.md`
3. If the request is downstream, immediately route conceptually to:
   - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-downstream-pipeline/SKILL.md`
4. If the user has just completed upstream work and now wants downstream analysis, reuse validated upstream artifacts instead of forcing the user to restate them.
5. Ask questions progressively. Do not front-load later-stage parameters.
6. Before running any substantial command, show the exact final command and wait for confirmation.
7. Execute confirmed commands directly in the shell or in the required `tmux` session.
8. Do not create temporary `.sh` script files for routine command execution unless the user explicitly asks for a reusable script.
9. Confirm one top-level results root once per workflow or continuation, then reuse it across later skills and steps.
10. After the top-level results root has been confirmed, derive step-specific output subdirectories automatically under that root instead of repeatedly asking the user to provide new output directories.
11. Show each resolved output path to the user for confirmation before execution.
12. Do not silently rerun long jobs if outputs may already exist.
13. If required inputs or dependencies are missing, stop and report clearly.
14. Some steps may take substantial time depending on data size, tool choice, environment, and parameters; a long runtime or a period without visible output does not necessarily mean failure.
15. For potentially long-running steps, prefer resumable execution patterns when possible, keep logs, and do not restart the same step unless failure is confirmed.

## Entry Procedure

At the start of a new request:

1. Read `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/README.md`.
2. Confirm the overall intent:
   - upstream only
   - downstream only
   - or full workflow / cross-stage continuation
3. Confirm one top-level results root under `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/results` and reuse it for later skill-specific subdirectories.
4. If the request starts from raw data, also confirm:
   - platform: `ONT` or `PacBio`
   - sample or run name
   - raw read file path
   - reference genome FASTA path
   - annotation GTF/GFF path
5. If the request starts from upstream outputs, confirm which validated outputs already exist, such as:
   - gene-by-cell matrix
   - transcript-by-cell or isoform-by-cell matrix
   - transcript annotation or extended annotation
   - Seurat RDS objects
   - cluster labels or reduced dimensions if they already exist
6. Validate that all declared input files exist before building commands.

Only after that should this skill hand off to the appropriate workflow entry skill.

## Routing Rules

### Route To Upstream

Route to the upstream entry skill if the user is asking for:

- ONT or PacBio preprocessing from raw reads
- BLAZE demultiplexing
- minimap2 alignment
- samtools BAM preparation
- IsoQuant, bambu, or FLAIR quantification
- transcript-by-cell or gene-by-cell matrix generation

Use:

- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-upstream-agent/SKILL.md`

### Route To Downstream

Route to the downstream entry skill if the user is asking for:

- resource setup for downstream analysis
- split gene count matrix by sample
- empty droplet removal or downstream QC
- sample merging or batch correction
- isoform assay integration into Seurat objects
- annotation
- isoforms of interest exploration
- isoform classification
- isoform diversity profiling
- differential transcript usage

Use:

- `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-downstream-pipeline/SKILL.md`

### Route Across Stages

If the user starts with upstream and then continues into downstream:

- first finish or validate the needed upstream outputs
- then explicitly state that the workflow is moving into downstream analysis
- then hand off to the downstream entry skill
- reuse upstream artifacts already produced in the current project context

## Project References And Working Paths

- project overview:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/README.md`
- upstream entry:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-upstream-agent/SKILL.md`
- downstream entry:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/.pi/skills/lr-scrna-downstream-pipeline/SKILL.md`
- input data directory:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/data`
- output root:
  - `/home4/zcz/zcz_proj1/IsoAgent/agent_proj/results`

## What To Deliver

For a normal successful run, this top-level skill should deliver:

- clear routing into the correct workflow layer
- reuse of existing validated project artifacts when possible
- the correct entry skill for the current stage
- stage-appropriate commands and outputs after handoff

## Response Style

Be operational and explicit:

- say whether the request is upstream, downstream, or cross-stage
- say which workflow entry skill is being used next
- say which validated inputs are being reused
- do not duplicate detailed module logic that already lives in the specialized entry skills
