#!/usr/bin/env python3
"""bj-expression-wrapper: run BioSkryb's BJ-Expression RNA-seq pipeline from ClawBio.

BJ-Expression takes already-demultiplexed FASTQ pairs (one pair per biosample)
and runs seqtk downsampling, cutadapt/fastp trimming, Salmon, STAR + HTSeq,
Qualimap and MultiQC.

This wrapper validates the input sheet, checks the tools, builds one fixed
`nextflow run` command, runs it, and summarises what the pipeline published.
It does NOT contain or download any BioSkryb code: the user supplies their own
checkout of the pinned tag, and must accept BioSkryb's terms to run it.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_DIR = Path(__file__).resolve().parent
DEMO_CSV = SKILL_DIR / "examples" / "demo_samplesheet.csv"

PROG = "bj-expression-wrapper"
PIPELINE_REPO = "https://github.com/BioSkryb/bj-expression.git"
DEFAULT_PIPELINE_VERSION = "2.0.0"
# Nextflow 26.x cannot parse bj-expression 2.0.0's config: it includes two
# profile files (batch_dev_gpu.config, batch_prod_gpu.config) that are not in
# the public repo. 25.04.6 and 24.10.5 both parse it (config parse only).
DEFAULT_NEXTFLOW_VERSION = "25.04.6"
LICENSE_URL = "https://www.bioskryb.com/basejumper/local-license/"
# Fixed so the published directory is predictable: <publish_dir>_<timestamp>.
RUN_TIMESTAMP = "clawbio"
GENOMES = ("GRCh38", "GRCm39")  # the two the BJ-Expression README documents
PROFILES = ("docker", "singularity")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
REMOTE_PREFIXES = ("s3://", "http://", "https://", "gs://", "az://", "ftp://")

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# Published by workflow/rnaseq.nf under <publish_dir>_<timestamp>/ (required
# ones must exist for a run to count as a success).
REQUIRED_OUTPUTS = {
    "multiqc_report": "multiqc/multiqc_report.html",
    "pipeline_metrics_summary": "secondary_analyses/secondary_metrics/pipeline_metrics_summary.csv",
    "salmon_gene_counts": "secondary_analyses/quantification_salmon/df_gene_counts_salmon.tsv",
    "htseq_gene_counts": "secondary_analyses/quantification_htseq/df_gene_counts_starhtseq.tsv",
}
OPTIONAL_OUTPUTS = {
    "pipeline_metrics_summary_percents": "secondary_analyses/secondary_metrics/pipeline_metrics_summary_percents.csv",
    "dynamic_range_expression": "secondary_analyses/secondary_metrics/df_dynamicrange_expression.tsv",
}
ALIGNMENT_DIR = "secondary_analyses/alignment_htseq"

REQUIRED_COLUMNS = ("biosampleName", "read1", "read2")


class WrapperError(Exception):
    """A user-facing failure; printed as one line and exits 1."""


# ──────────────────────────────────────────────
# Input sheet
# ──────────────────────────────────────────────

def is_safe_sample_name(name: str) -> bool:
    """Sample names become file names and Nextflow tags: keep them plain."""
    return bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name)) and ".." not in name


def _resolve_fastq(value: str, base_dir: Path) -> tuple[str, bool]:
    """Return (path-or-url, is_remote); local paths become absolute."""
    if value.lower().startswith(REMOTE_PREFIXES):
        return value, True
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = base_dir / p
    return str(p.resolve()), False


def validate_input_csv(path: Path, require_groups: bool,
                       base_dir: Path | None = None) -> list[dict]:
    """Parse and validate the sample sheet; return normalised rows.

    Relative FASTQ paths resolve against `base_dir` (default: the current
    directory, as Nextflow does). Remote URLs are accepted without a check.
    """
    path = Path(path)
    if not path.is_file():
        raise WrapperError(f"input sheet not found: {path}")
    base = Path(base_dir) if base_dir else Path.cwd()

    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        header = [h.strip() for h in (reader.fieldnames or [])]
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise WrapperError(
                f"input sheet is missing required column(s): {', '.join(missing)} "
                f"(need biosampleName,read1,read2[,groups])")
        raw = [{(k or "").strip(): (v or "").strip() for k, v in r.items()}
               for r in reader]
    if not raw:
        raise WrapperError("input sheet has no samples (header only)")

    rows, seen = [], set()
    for i, r in enumerate(raw, start=2):  # line numbers, header is line 1
        name = r["biosampleName"]
        if not is_safe_sample_name(name):
            raise WrapperError(
                f"line {i}: biosampleName {name!r} is not allowed; use letters, "
                "digits, '.', '_' and '-' only, starting with a letter or digit")
        if name in seen:
            raise WrapperError(f"line {i}: duplicate biosampleName {name!r}")
        seen.add(name)

        fq, remote = [], False
        for col in ("read1", "read2"):
            value = r[col]
            if not value:
                raise WrapperError(f"line {i}: {col} is empty for {name}")
            if not value.lower().endswith(FASTQ_SUFFIXES):
                raise WrapperError(
                    f"line {i}: {col} {value!r} is not a FASTQ "
                    f"({', '.join(FASTQ_SUFFIXES)})")
            resolved, is_remote = _resolve_fastq(value, base)
            if not is_remote and not Path(resolved).is_file():
                raise WrapperError(f"line {i}: {col} file not found: {resolved}")
            remote = remote or is_remote
            fq.append(resolved)
        if fq[0] == fq[1]:
            raise WrapperError(
                f"line {i}: read1 and read2 are the same file for {name}")

        rows.append({"biosampleName": name, "read1": fq[0], "read2": fq[1],
                     "groups": r.get("groups", ""), "_remote": remote})

    if require_groups and any(not r["groups"] for r in rows):
        raise WrapperError(
            "the 'groups' column is required for every sample when QC plots are "
            "on (per the BJ-Expression README). Add it, or pass --skip-qc-plots.")
    return rows


def write_normalised_csv(rows: list[dict], dest: Path) -> Path:
    """Write the sheet the pipeline will read: absolute paths, known columns."""
    with_groups = any(r["groups"] for r in rows)
    cols = [*REQUIRED_COLUMNS, *(["groups"] if with_groups else [])]
    with open(dest, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])
    return dest


# ──────────────────────────────────────────────
# Pipeline source, tools
# ──────────────────────────────────────────────

def _git(path: Path, *args: str) -> str:
    try:
        proc = subprocess.run(["git", *args], cwd=str(path), capture_output=True,
                              text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise WrapperError(f"could not run git in {path}: {exc}") from exc
    return proc.stdout.strip() if proc.returncode == 0 else ""


def check_pipeline_source(pipeline_dir: Path | None, version: str,
                          allow_dirty: bool) -> dict:
    """The checkout must be the pinned tag, unmodified (unless allowed)."""
    how_to = (
        f"Get the pipeline yourself, from BioSkryb, at the pinned tag:\n"
        f"  git clone --branch {version} --depth 1 {PIPELINE_REPO}\n"
        f"then pass --pipeline-dir <that folder>. Using it is subject to "
        f"BioSkryb's terms: {LICENSE_URL}")
    if pipeline_dir is None:
        raise WrapperError(f"--pipeline-dir is required.\n{how_to}")
    d = Path(pipeline_dir).expanduser().resolve()
    if not (d / "main.nf").is_file() or not (d / "nextflow.config").is_file():
        raise WrapperError(f"no BJ-Expression checkout at {d} (main.nf or "
                           f"nextflow.config missing).\n{how_to}")
    tag = _git(d, "describe", "--tags", "--exact-match")
    commit = _git(d, "rev-parse", "HEAD")
    if tag != version:
        raise WrapperError(
            f"{d} is at {tag or 'an untagged commit'}, not the pinned {version}. "
            f"Run: git -C {d} checkout {version}")
    dirty = bool(_git(d, "status", "--porcelain"))
    if dirty and not allow_dirty:
        raise WrapperError(
            f"{d} has uncommitted changes, so the run would not match tag "
            f"{version}. Revert them, or pass --allow-dirty-pipeline "
            "(recorded in the result).")
    return {"repo": PIPELINE_REPO, "path": str(d), "version": version,
            "commit": commit, "dirty": dirty, "license_url": LICENSE_URL}


def nextflow_env(nxf_version: str) -> dict:
    """Environment for every nextflow call: pin NXF_VER unless told 'system'."""
    env = dict(os.environ)
    if nxf_version != "system":
        env["NXF_VER"] = nxf_version
    return env


def _version_line(cmd: list[str], pattern: str, env: dict | None = None) -> str:
    try:
        # the first call with a pinned NXF_VER downloads that Nextflow release
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                              env=env)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    m = re.search(pattern, proc.stdout + proc.stderr)
    return m.group(1) if m else "unknown"


def preflight_tools(profile: str, nxf_version: str = DEFAULT_NEXTFLOW_VERSION) -> dict:
    """Check Java, Nextflow and the container backend; report every problem."""
    problems, tools = [], {}
    if shutil.which("java"):
        tools["java"] = _version_line(["java", "-version"], r'version "([^"]+)"')
    else:
        problems.append(
            "java not found on PATH. Nextflow needs a JDK "
            "(for example: conda install -c conda-forge openjdk).")
    if shutil.which("nextflow"):
        tools["nextflow"] = _version_line(["nextflow", "-version"],
                                          r"version\s+(\S+)",
                                          env=nextflow_env(nxf_version))
    else:
        problems.append("nextflow not found on PATH. Install it from "
                        "https://www.nextflow.io/docs/latest/install.html")

    if profile == "docker":
        if not shutil.which("docker"):
            problems.append(
                "docker not found on PATH. Install Docker (needs a Carolina "
                "Cloud VM, not a container), or use --profile singularity.")
        else:
            ok = subprocess.run(["docker", "info"], capture_output=True).returncode == 0
            if ok:
                tools["docker"] = "available"
            else:
                problems.append("docker is installed but its daemon is not "
                                "reachable (`docker info` failed).")
    else:
        sif = shutil.which("apptainer") or shutil.which("singularity")
        if sif:
            tools["singularity"] = Path(sif).name
        else:
            problems.append("--profile singularity needs apptainer or "
                            "singularity on PATH; neither was found.")
    if problems:
        raise WrapperError("preflight failed:\n  - " + "\n  - ".join(problems))
    return tools


# ──────────────────────────────────────────────
# Params and command
# ──────────────────────────────────────────────

def build_params(args: argparse.Namespace, csv_path: Path, out_dir: Path) -> dict:
    """Required values plus only the overrides the user asked for.

    Anything not listed keeps the pipeline's own default; passing our own
    defaults would silently change BioSkryb's method.
    """
    params: dict = {
        "input_csv": str(csv_path),
        "publish_dir": str(out_dir / "results"),
        "timestamp": RUN_TIMESTAMP,
        "genome": args.genome,
    }
    for attr, key in (("n_reads", "n_reads"), ("read_length", "read_length"),
                      ("min_reads", "min_reads"), ("seed", "seqtk_sample_seed"),
                      ("instrument", "instrument"), ("max_cpus", "max_cpus"),
                      ("max_memory", "max_memory"), ("genomes_base", "genomes_base")):
        value = getattr(args, attr)
        if value is not None:
            params[key] = value
    for flag in ("skip_subsampling", "skip_cutadapt", "skip_qc_plots"):
        if getattr(args, flag):
            params[flag] = True
    return params


def build_command(pipeline_dir: str, params_file: Path, work_dir: Path,
                  profile: str, resume: bool,
                  extra_config: Path | None = None,
                  trace: Path | None = None) -> list[str]:
    cmd = ["nextflow", "run", str(Path(pipeline_dir) / "main.nf"),
           "-params-file", str(params_file), "-work-dir", str(work_dir),
           "-ansi-log", "false"]
    if trace is not None:
        cmd += ["-with-trace", str(trace)]
    if extra_config is not None:
        cmd += ["-c", str(extra_config)]
    if profile != "docker":  # docker is enabled by the pipeline's own config
        cmd += ["-profile", profile]
    if resume:
        cmd.append("-resume")
    return cmd


# ──────────────────────────────────────────────
# Execution and outputs
# ──────────────────────────────────────────────

def run_nextflow(cmd: list[str], run_dir: Path,
                 env: dict | None = None) -> tuple[int, str]:
    """Run in run_dir (Nextflow litters its cwd); logs go to files, tail returned."""
    run_dir.mkdir(parents=True, exist_ok=True)
    out_log, err_log = run_dir / "nextflow.stdout.log", run_dir / "nextflow.stderr.log"
    with open(out_log, "w") as so, open(err_log, "w") as se:
        proc = subprocess.run(cmd, cwd=run_dir, stdout=so, stderr=se, env=env)
    tail = "\n".join(err_log.read_text(errors="replace").splitlines()[-40:])
    return proc.returncode, tail


FAILED_STATES = ("FAILED", "ABORTED")


def parse_trace(path: Path) -> dict | None:
    """Summarise a Nextflow trace by each task's FINAL state.

    BJ-Expression sets errorStrategy to retry twice and then *ignore*, so a run
    can exit 0 with steps missing. A task name whose last row is FAILED/ABORTED
    was ignored; one that failed and later completed was only retried.
    Returns None if the trace is missing or has no name/status columns.
    """
    path = Path(path)
    if not path.is_file():
        return None
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"name", "status"} <= set(reader.fieldnames):
            return None
        history: dict[str, list[str]] = {}
        for row in reader:
            if row.get("name"):
                history.setdefault(row["name"], []).append((row.get("status") or "").strip())
    ignored = sorted(n for n, h in history.items() if h[-1] in FAILED_STATES)
    retried = sorted(n for n, h in history.items()
                     if h[-1] not in FAILED_STATES and any(x in FAILED_STATES for x in h))
    return {"tasks": len(history), "ignored_failures": ignored, "retried_tasks": retried}


def collect_outputs(results_dir: Path, samples: list[str]) -> dict:
    found, missing = {"results_dir": str(results_dir)}, []
    for key, rel in REQUIRED_OUTPUTS.items():
        if (results_dir / rel).is_file():
            found[key] = str(results_dir / rel)
        else:
            missing.append(rel)
    for key, rel in OPTIONAL_OUTPUTS.items():
        if (results_dir / rel).is_file():
            found[key] = str(results_dir / rel)
    bams = sorted(p.name[:-4] for p in (results_dir / ALIGNMENT_DIR).glob("*.bam")) \
        if (results_dir / ALIGNMENT_DIR).is_dir() else []
    found["aligned_samples"] = bams
    return {"outputs": found, "missing_outputs": missing,
            "samples_missing_from_outputs": [s for s in samples if s not in bams]}


# ──────────────────────────────────────────────
# Report and bundle
# ──────────────────────────────────────────────

def write_report(out_dir: Path, result: dict) -> Path:
    status = result["status"]
    lines = [
        "# bj-expression-wrapper Report", "",
        f"**Date**: {result['generated_at']}",
        f"**Status**: {status}",
        f"**Genome**: {result['params']['genome']}",
        f"**Samples**: {len(result['samples'])} ({', '.join(result['samples'])})",
    ]
    pipe = result.get("pipeline")
    if pipe:
        lines.append(f"**Pipeline**: BioSkryb bj-expression {pipe['version']} "
                     f"(commit {pipe['commit'][:12] or 'unknown'}"
                     f"{', modified' if pipe['dirty'] else ''})")
    if result.get("tools"):
        lines.append("**Tools**: " + ", ".join(f"{k} {v}" for k, v in result["tools"].items()))
    lines.append("")

    if status == "demo_plan_only":
        lines += ["> **Plan only.** This demo validated a synthetic sample sheet and "
                  "wrote the plan. Nothing was run: no pipeline, Nextflow or "
                  "container was used.", ""]
    elif status == "check_passed":
        lines += ["> **Check only.** The inputs, pipeline checkout and tools passed; "
                  "the pipeline was not run.", ""]

    outs = result.get("outputs") or {}
    if status in ("success", "incomplete") and outs:
        lines += ["## Outputs", "", "| Output | Path |", "|---|---|"]
        for key in (*REQUIRED_OUTPUTS, *OPTIONAL_OUTPUTS):
            if key in outs:
                lines.append(f"| {key} | `{outs[key]}` |")
        lines.append("")
    if result.get("missing_outputs"):
        lines += ["**Expected but not produced:** "
                  + ", ".join(f"`{m}`" for m in result["missing_outputs"]), ""]
    dropped = result.get("samples_missing_from_outputs") or []
    if dropped and status in ("success", "incomplete"):
        lines += [f"> **Warning:** no alignment was published for: {', '.join(dropped)}. "
                  "BJ-Expression filters out samples with fewer reads than `min_reads` "
                  "before quantification, so they will be missing from the count "
                  "tables. Check the pipeline log before assuming a sample was analysed.",
                  ""]
    if status == "failed":
        lines += ["## Nextflow stderr (tail)", "", "```", result["stderr_tail"], "```", ""]

    ignored_tasks = result.get("ignored_failures")
    if ignored_tasks:
        lines += ["## Tasks that failed and were ignored", "",
                  "BJ-Expression retries a failed task twice and then ignores it, so the "
                  "run continues and can exit 0. These tasks did not finish:", ""]
        lines += [f"- `{t}`" for t in ignored_tasks]
        lines += ["", "Anything downstream of them may be missing or partial. Check "
                  "`nextflow_run/trace.tsv` and the task work directories.", ""]
    elif ignored_tasks is None and status in ("success", "success_with_ignored_failures",
                                              "incomplete"):
        lines += ["> **No Nextflow trace was available**, so ignored task failures could "
                  "not be checked. Do not read this status as proof every step finished.", ""]

    overridden = sorted(set(result["params"]) - {"input_csv", "publish_dir",
                                                 "timestamp", "genome"})
    lines += ["## Parameters", "",
              "Only these were set by the wrapper beyond the required ones: "
              + (", ".join(overridden) if overridden else "none")
              + ". Everything else is BioSkryb's default for the pinned tag.", "",
              "## Pipeline licence", "",
              f"BJ-Expression is BioSkryb software, used under BioSkryb's terms "
              f"({LICENSE_URL}). This wrapper contains none of its code.", "",
              f"*{DISCLAIMER}*", ""]
    path = out_dir / "report.md"
    path.write_text("\n".join(lines))
    return path


def write_bundle(out_dir: Path, command: str | None,
                 extra_config: Path | None = None) -> None:
    write_commands_sh(out_dir, command or "# plan only: no command was run")
    write_environment_yml(
        out_dir, "clawbio-bj-expression-wrapper", pip_deps=[],
        conda_deps=["openjdk", "nextflow"],
        channels=["conda-forge", "bioconda", "nodefaults"])
    repro = out_dir / "reproducibility"
    write_checksums(
        [out_dir / "report.md", out_dir / "result.json", repro / "commands.sh",
         repro / "params.json", repro / "samplesheet.normalised.csv",
         *([extra_config] if extra_config else [])],
        out_dir, anchor=out_dir)


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

def _int(minimum: int, what: str):
    def conv(value: str) -> int:
        try:
            n = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{what} must be an integer")
        if n < minimum:
            raise argparse.ArgumentTypeError(f"{what} must be >= {minimum}")
        return n
    return conv


def _nxf_version(value: str) -> str:
    if value == "system" or re.fullmatch(r"\d+\.\d+\.\d+", value):
        return value
    raise argparse.ArgumentTypeError(
        "--nextflow-version must look like 25.04.6 (or 'system' to use the "
        "installed Nextflow unchanged)")


def _memory(value: str) -> str:
    if not re.fullmatch(r"\d+(?:\.\d+)?[. ]?(?:KB|MB|GB|TB)", value):
        raise argparse.ArgumentTypeError("--max-memory looks like 50.GB or 512.MB")
    return value


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=PROG,
        description="Run BioSkryb BJ-Expression (RNA-seq from demultiplexed FASTQ pairs).")
    p.add_argument("--input", type=Path,
                   help="CSV: biosampleName,read1,read2[,groups]; relative paths "
                        "resolve against the current directory")
    p.add_argument("--output", type=Path, help="Output directory")
    p.add_argument("--demo", action="store_true",
                   help="Validate a bundled synthetic sheet and write the plan; runs nothing")
    p.add_argument("--check", action="store_true",
                   help="Validate inputs, pipeline checkout and tools; do not run")
    p.add_argument("--pipeline-dir", type=Path,
                   help="Your checkout of BioSkryb/bj-expression at the pinned tag")
    p.add_argument("--pipeline-version", default=DEFAULT_PIPELINE_VERSION,
                   help=f"Tag the checkout must be at (default {DEFAULT_PIPELINE_VERSION})")
    p.add_argument("--accept-bioskryb-terms", action="store_true",
                   help=f"Confirm you accept BioSkryb's terms ({LICENSE_URL}); needed to run")
    p.add_argument("--allow-dirty-pipeline", action="store_true",
                   help="Allow uncommitted changes in the checkout (recorded)")
    p.add_argument("--genome", choices=GENOMES, default="GRCh38")
    p.add_argument("--profile", choices=PROFILES, default="docker",
                   help="Container backend (default docker)")
    p.add_argument("--n-reads", type=_int(1, "--n-reads"),
                   help="Reads to sample per sample (pipeline default if omitted)")
    p.add_argument("--read-length", type=_int(1, "--read-length"),
                   help="Read length to trim to (pipeline default if omitted)")
    p.add_argument("--min-reads", type=_int(0, "--min-reads"),
                   help="Minimum reads for a sample to be analysed (pipeline default if omitted)")
    p.add_argument("--seed", type=_int(0, "--seed"),
                   help="seqtk sample seed inside the pipeline (pipeline default if omitted)")
    p.add_argument("--skip-subsampling", action="store_true")
    p.add_argument("--skip-cutadapt", action="store_true")
    p.add_argument("--skip-qc-plots", action="store_true")
    p.add_argument("--instrument", help="Sequencer, e.g. NovaSeq (selects two-colour chemistry)")
    p.add_argument("--max-cpus", type=_int(1, "--max-cpus"))
    p.add_argument("--max-memory", type=_memory)
    p.add_argument("--genomes-base", help="Mirror of the reference bundle (default: BioSkryb's public S3)")
    p.add_argument("--resume", action="store_true", help="Pass -resume to Nextflow")
    p.add_argument("--strict", action="store_true",
                   help="Exit 1 if any task failed and was ignored by the pipeline")
    p.add_argument("--nextflow-version", type=_nxf_version,
                   default=DEFAULT_NEXTFLOW_VERSION,
                   help=f"Nextflow release to run (default {DEFAULT_NEXTFLOW_VERSION}; "
                        "'system' = whatever is installed). Newer 26.x releases "
                        "cannot parse this pipeline's config.")
    p.add_argument("--extra-config", type=Path,
                   help="Extra Nextflow config passed with -c (for example container "
                        "settings for hosts without Docker); copied into the bundle")
    return p


def run(args: argparse.Namespace) -> int:
    if args.output is None:
        raise WrapperError("--output is required")
    if args.demo:
        mode = "demo"
    elif args.check:
        mode = "check"
    else:
        mode = "run"
        if args.input is None:
            raise WrapperError("--input is required (or use --demo / --check)")
        if not args.accept_bioskryb_terms:
            raise WrapperError(
                "BJ-Expression is BioSkryb software under BioSkryb's terms "
                f"({LICENSE_URL}). Read them, then re-run with "
                "--accept-bioskryb-terms to confirm.")
    if mode == "check" and args.input is None:
        raise WrapperError("--input is required with --check")

    # Nextflow runs from <output>/nextflow_run, so every path must be absolute.
    out_dir: Path = args.output.expanduser().resolve()
    if (out_dir / "result.json").exists() and not (mode == "run" and args.resume):
        raise WrapperError(
            f"{out_dir} already holds results; choose a new --output "
            "(or pass --resume to continue a Nextflow run).")

    if mode == "demo":
        csv_in, base = DEMO_CSV, DEMO_CSV.parent
    else:
        csv_in, base = args.input, None
    rows = validate_input_csv(csv_in, require_groups=not args.skip_qc_plots,
                              base_dir=base)
    samples = [r["biosampleName"] for r in rows]

    pipeline, tools = None, {}
    if mode != "demo":
        pipeline = check_pipeline_source(args.pipeline_dir, args.pipeline_version,
                                         args.allow_dirty_pipeline)
        tools = preflight_tools(args.profile, args.nextflow_version)

    if args.extra_config is not None and not args.extra_config.is_file():
        raise WrapperError(f"--extra-config file not found: {args.extra_config}")

    out_dir.mkdir(parents=True, exist_ok=True)
    repro = out_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    kept_config = None
    if args.extra_config is not None:
        kept_config = repro / "extra.config"
        shutil.copyfile(args.extra_config, kept_config)
    csv_path = write_normalised_csv(rows, repro / "samplesheet.normalised.csv")
    params = build_params(args, csv_path, out_dir)
    (repro / "params.json").write_text(json.dumps(params, indent=2) + "\n")

    run_dir = out_dir / "nextflow_run"
    cmd = build_command(pipeline["path"] if pipeline else "<path-to-bj-expression-checkout>",
                        repro / "params.json", out_dir / "work", args.profile, args.resume,
                        kept_config, run_dir / "trace.tsv")
    pin = "" if args.nextflow_version == "system" else f"NXF_VER={args.nextflow_version} "
    command_text = f"cd {shlex.quote(str(run_dir))} && {pin}{shlex.join(cmd)}"

    result = {
        "skill": PROG,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": {"demo": "demo_plan_only", "check": "check_passed", "run": "running"}[mode],
        "terms_accepted": bool(args.accept_bioskryb_terms),
        "pipeline": pipeline,
        "tools": tools,
        "profile": args.profile,
        "nextflow_version_requested": args.nextflow_version,
        "extra_config": str(kept_config) if kept_config else None,
        "params": params,
        "samples": samples,
        "command": command_text,
        "disclaimer": DISCLAIMER,
    }
    exit_code = 0
    if mode == "run":
        code, tail = run_nextflow(cmd, run_dir, nextflow_env(args.nextflow_version))
        result["nextflow_exit_code"] = code
        result["stderr_tail"] = tail
        found = collect_outputs(Path(params["publish_dir"] + "_" + RUN_TIMESTAMP), samples)
        result.update(found)
        trace = parse_trace(run_dir / "trace.tsv")
        result["trace"] = ({"available": True, "tasks": trace["tasks"],
                            "retried_tasks": trace["retried_tasks"]} if trace
                           else {"available": False, "tasks": 0, "retried_tasks": []})
        result["ignored_failures"] = trace["ignored_failures"] if trace else None
        if code != 0:
            result["status"], exit_code = "failed", 1
        elif found["missing_outputs"]:
            result["status"], exit_code = "incomplete", 1
        elif result["ignored_failures"]:
            result["status"] = "success_with_ignored_failures"
            exit_code = 1 if args.strict else 0
        else:
            result["status"] = "success"

    (out_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    write_report(out_dir, result)
    write_bundle(out_dir, None if mode == "demo" else command_text, kept_config)

    if result["status"] == "failed":
        print(f"{PROG}: error: Nextflow exited {result['nextflow_exit_code']}\n"
              f"{result['stderr_tail']}", file=sys.stderr)
    elif result["status"] == "incomplete":
        print(f"{PROG}: error: Nextflow exited 0 but did not publish: "
              + ", ".join(result["missing_outputs"]), file=sys.stderr)
    else:
        print(f"{PROG}: {result['status']} -> {out_dir / 'report.md'}")
    if result["status"] == "success_with_ignored_failures":
        names = ", ".join(result["ignored_failures"])
        msg = f"tasks failed and were ignored by the pipeline: {names}"
        if exit_code:
            print(f"{PROG}: error: {msg} (--strict)", file=sys.stderr)
        else:
            print(f"{PROG}: WARNING: {msg}")
    return exit_code


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except WrapperError as exc:
        print(f"{PROG}: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
