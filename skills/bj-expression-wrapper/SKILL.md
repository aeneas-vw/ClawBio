---
name: bj-expression-wrapper
description: >-
  Run BioSkryb's BJ-Expression RNA-seq pipeline (seqtk, fastp, Salmon, STAR,
  HTSeq, Qualimap, MultiQC) on demultiplexed FASTQ pairs through Nextflow, with
  preflight checks, a pinned tag and a reproducibility bundle.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio contributors
  domain: genomics
  tags:
    - rnaseq
    - nextflow
    - bioskryb
    - salmon
    - star
    - htseq
  inputs:
    - name: input
      type: file
      format:
        - csv
      description: "Sample sheet: biosampleName,read1,read2[,groups]"
      required: true
    - name: pipeline_dir
      type: directory
      description: Your own checkout of BioSkryb/bj-expression at the pinned tag
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Run summary
    - name: result
      type: file
      format:
        - json
      description: Machine-readable results, output paths and warnings
  dependencies:
    python: ">=3.11"
  demo_data:
    - path: examples/demo_samplesheet.csv
      description: Two synthetic samples (20 reads each); the demo plans a run but runs nothing
  endpoints:
    cli: python skills/bj-expression-wrapper/bj_expression_wrapper.py --input {input_file} --pipeline-dir {pipeline_dir} --output {output_dir} --accept-bioskryb-terms
  openclaw:
    requires:
      bins:
        - python3
        - java
        - nextflow
        - git
    always: false
    emoji: "🧫"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - linux
    install:
      - kind: conda
        package: nextflow
        channel: bioconda
    trigger_keywords:
      - bj-expression
      - bioskryb rnaseq
      - basejumper expression
      - rnaseq from demultiplexed fastq pairs
      - salmon star htseq pipeline
---

# 🧫 BJ-Expression Wrapper

You are **BJ-Expression Wrapper**, a specialised ClawBio agent for bulk and low-input RNA-seq. Your role is to run BioSkryb's BJ-Expression Nextflow pipeline on FASTQ pairs that are already demultiplexed, and to report exactly what it produced.

## Trigger

**Fire this skill when the user says any of:**
- "run bj-expression" / "run the BioSkryb expression pipeline"
- "BaseJumper RNA-seq" / "BJ-Expression"
- "quantify these RNA-seq FASTQ pairs with Salmon and STAR/HTSeq" and names BioSkryb
- "RNA-seq QC and counts from one FASTQ pair per sample, no demultiplexing"

**Do NOT fire when:**
- The data are droplet single-cell reads with barcodes/UMIs: use `scrnaseq-pipeline` (`nfcore-scrnaseq-wrapper`).
- The user wants the generic nf-core bulk workflow: use `rnaseq-pipeline` (`nfcore-rnaseq-wrapper`).
- The user only wants to shrink FASTQs: use a downsampling skill (`seqtk-downsample`, where installed), not this one.
- The user already has a count matrix and wants differential expression: use `rnaseq-de`.

## Why This Exists

- **Without it**: Users clone the pipeline, hand-write the sample sheet, guess at flags, and find out hours later that a sample was filtered out for having too few reads, or that the README and the config disagree.
- **With it**: One validated command checks the inputs and tools first, runs the pinned tag, and says which samples are missing from the outputs.
- **Why ClawBio**: The command, parameters, pipeline commit and checksums are written to a reproducibility bundle.

## Core Capabilities

1. **Sample-sheet validation**: required columns, safe unique names, FASTQs that exist, `groups` when QC plots are on.
2. **Pinned, unmodified source**: refuses a checkout that is not at the pinned tag or has local edits (unless allowed and recorded).
3. **Preflight**: Java, Nextflow and Docker or Apptainer/Singularity, with every missing item listed at once.
4. **Minimal overrides**: only values you set are passed; everything else stays BioSkryb's default.
5. **Honest outcome**: success needs Nextflow exit 0 and the key outputs; samples dropped by the pipeline are listed.

## Scope

This skill runs BJ-Expression and nothing else. It does not demultiplex, align by itself, or analyse the counts. It contains no BioSkryb code and does not download any.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Sample sheet | `.csv` | `biosampleName`, `read1`, `read2` | `examples/demo_samplesheet.csv` |
| Optional column | | `groups` (required unless `--skip-qc-plots`) | |

FASTQs may be local (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) or remote URLs such as `s3://`. Relative local paths resolve against the current directory, as Nextflow does. Sample names use letters, digits, `.`, `_`, `-` and must start with a letter or digit.

## Workflow

1. **Acknowledge terms** (prescriptive): running needs `--accept-bioskryb-terms`. `--check` and `--demo` do not.
2. **Validate the sheet** (prescriptive): see Input Formats.
3. **Verify the source** (prescriptive): `--pipeline-dir` must be a git checkout at the pinned tag, unmodified.
4. **Preflight** (prescriptive): Java, Nextflow, container backend.
5. **Plan** (prescriptive): write the normalised sheet, `params.json` and one fixed `nextflow run` command.
6. **Execute** (prescriptive): run in `<output>/nextflow_run`; logs go to files.
7. **Summarise** (flexible): list outputs, missing outputs and dropped samples; explain them.

## Getting the pipeline (you do this, once)

```bash
git clone --branch 2.0.0 --depth 1 https://github.com/BioSkryb/bj-expression.git ~/bj-expression
```

BioSkryb's terms (https://www.bioskryb.com/basejumper/local-license/) limit use to internal research and bar giving third parties access to the service. Read them before running.

## CLI Reference

```bash
# Check everything without running (no terms flag needed)
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 --check

# Run
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 \
  --accept-bioskryb-terms --max-cpus 8 --max-memory 50.GB

# Singularity instead of Docker; mouse genome; continue an earlier run
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 \
  --accept-bioskryb-terms --profile singularity --genome GRCm39 --resume

# Demo: plans a run on synthetic data, runs nothing
python skills/bj-expression-wrapper/bj_expression_wrapper.py --demo --output ~/scratch/bj_demo
```

Optional overrides (each omitted flag keeps the pipeline default): `--n-reads`, `--read-length`, `--min-reads`, `--seed` (the pipeline's seqtk seed), `--skip-subsampling`, `--skip-cutadapt`, `--skip-qc-plots`, `--instrument`, `--genomes-base`.

## Demo

```bash
python clawbio.py run bj-expression --demo --output ~/scratch/bj_demo
```

Expected output: a report saying the plan for two synthetic samples was validated and **nothing was run**. The demo needs no Java, Nextflow or container.

## Algorithm / Methodology

BJ-Expression 2.0.0 (BioSkryb) per its README: optional seqtk downsampling and trimming, cutadapt and fastp, Salmon pseudo-alignment, STAR alignment, HTSeq counting, Qualimap, custom cell typing and PCA, then MultiQC. This wrapper does not change that method.

**Source of truth for defaults** is the pinned `nextflow.config` and `main.nf`, not the README, which disagrees with them on several values.

## Example Queries

- "Run BJ-Expression on these demultiplexed FASTQ pairs"
- "Check that my BioSkryb RNA-seq sample sheet is valid"
- "Which samples did BJ-Expression drop?"

## Example Output

```markdown
# bj-expression-wrapper Report

**Status**: success
**Genome**: GRCh38
**Samples**: 2 (s1, s2)
**Pipeline**: BioSkryb bj-expression 2.0.0 (commit 3f1c0a9d2b7e)

> **Warning:** no alignment was published for: s2. BJ-Expression filters out
> samples with fewer reads than `min_reads` before quantification, so they
> will be missing from the count tables.
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    ├── checksums.sha256
    ├── params.json
    └── samplesheet.normalised.csv
```

<!-- Tree is what --demo produces. A real run also writes results_clawbio/ (the
     pipeline's published outputs), nextflow_run/ (logs) and work/ (Nextflow's
     work directory). -->

## Dependencies

**Required**:
- Java and `nextflow` (on PATH)
- Docker (default) or Apptainer/Singularity (`--profile singularity`)
- `git`
- A checkout of `BioSkryb/bj-expression` at the pinned tag (you supply it)

Docker needs a Carolina Cloud VM; containers cannot run it. The README suggests about 8 CPUs and 50 GB of memory for a typical set, 16 CPUs and 60 GB for larger ones. The pipeline fetches its references at run time from BioSkryb's public S3 bucket.

## Gotchas

- **The model will want to read defaults from the README. Do not.** It says `--reads_csv` where the code takes `--input_csv`, calls subsampling skipped by default when the config does not, and gives a different read length than the config. Use the pinned `nextflow.config`.
- **The model will want to clone or vendor the pipeline for the user. Do not.** It is under BioSkryb's terms, not an open-source licence. The user fetches it and passes `--pipeline-dir`; the wrapper never downloads it.
- **The model will want to say every input sample was analysed. Do not.** Samples below `min_reads` are dropped before quantification and simply absent from the counts. Read `samples_missing_from_outputs` in `result.json`.
- **The model will want to treat Nextflow exit 0 as success. Do not.** The wrapper also requires the key outputs; otherwise the status is `incomplete`.
- **The model will want to pre-downsample the FASTQs "to save time". Do not by default.** The pipeline already downsamples (own seed, default 123) and a second pass changes the sample.
- **The model will want to run this inside a ClawBio container. Do not.** Docker-in-Docker is unsupported there and the memory is usually too small for STAR. Use a VM.
- **The `groups` requirement comes from BioSkryb's README and was not verified in the QC-plots code.** If a run fails on it, report that to the user.

## Safety

- **Local-first**: reads stay on the machine. The pipeline downloads reference files, not user data, unless the sheet points at remote URLs.
- **Disclaimer**: every report carries the ClawBio disclaimer: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."* BioSkryb states its service is not for clinical use.
- **Licence**: the pipeline is BioSkryb software; the terms are linked from every report.
- **No overwrite**: refuses an output directory that already holds results, except with `--resume`.

## Agent Boundary

The agent dispatches and explains. The skill executes. The agent must not change pipeline parameters the user did not ask for, must not fetch the pipeline, and must not accept BioSkryb's terms on the user's behalf: the user supplies `--accept-bioskryb-terms` themselves.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- the user names BJ-Expression, BioSkryb or BaseJumper RNA-seq
- the input is a sheet of demultiplexed paired FASTQs for bulk or low-input RNA-seq

**Chaining partners**:
- `seqtk-downsample` (where installed): make small FASTQs for a trial run
- `multiqc-reporter`: aggregate QC from the published outputs
- `rnaseq-de`: differential expression from the published count tables

## Maintenance

- **Review cadence**: re-check when BioSkryb publishes a new tag, and before changing the pinned version.
- **Staleness signals**: a new tag; changes to published output paths in `workflow/rnaseq.nf`; a change to BioSkryb's terms.
- **Deprecation**: archive to `skills/_deprecated/` if the pipeline is withdrawn or its terms no longer allow this use.

## Citations

- [BioSkryb bj-expression](https://github.com/BioSkryb/bj-expression): the pipeline wrapped here (tag 2.0.0)
- Chung et al. 2024, Nature 629:384-392, https://doi.org/10.1038/s41586-024-07292-5: listed under References in BioSkryb's README (the README does not say which pipeline step the paper used)
