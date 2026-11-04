"""Tests for the bj-expression-wrapper skill (red/green TDD).

BioSkryb's bj-expression pipeline needs Java, Nextflow and Docker/Singularity,
none of which CI has. These tests therefore use a fake `nextflow` (plus fake
`java` and `docker`) on a restricted PATH. The fake writes the output tree that
`workflow/rnaseq.nf` publishes. Nothing here runs the real pipeline.
"""

import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "bj_expression_wrapper.py"
LICENSE_URL = "https://www.bioskryb.com/basejumper/local-license/"
PIPELINE_TAG = "2.0.0"


def _load_module():
    spec = importlib.util.spec_from_file_location("bj_expression_wrapper", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bj_expression_wrapper"] = mod
    spec.loader.exec_module(mod)
    return mod


# ──────────────────────────────────────────────
# Fixtures: fake tools, fake pipeline checkout, fake FASTQs
# ──────────────────────────────────────────────

FAKE_NEXTFLOW = f"""#!{sys.executable}
import json, os, sys
args = sys.argv[1:]
log = os.environ.get("FAKE_NF_LOG")
if log:
    with open(log, "a") as fh:
        fh.write(json.dumps(args) + "\\n")
if args[:1] == ["-version"]:
    print("nextflow version 25.04.0 build 5000")
    sys.exit(0)
mode = os.environ.get("FAKE_NF_MODE", "ok")
if mode == "fail":
    print("boom: process STARALIGN failed", file=sys.stderr)
    sys.exit(3)
if mode == "empty":
    sys.exit(0)
params = json.load(open(args[args.index("-params-file") + 1]))
res = params["publish_dir"] + "_" + params["timestamp"]
sa = res + "/secondary_analyses"
drop = os.environ.get("FAKE_NF_DROP", "")
import csv
samples = [r["biosampleName"] for r in csv.DictReader(open(params["input_csv"]))]
files = [
    res + "/multiqc/multiqc_report.html",
    sa + "/secondary_metrics/pipeline_metrics_summary.csv",
    sa + "/secondary_metrics/pipeline_metrics_summary_percents.csv",
    sa + "/secondary_metrics/df_dynamicrange_expression.tsv",
    sa + "/quantification_salmon/df_gene_counts_salmon.tsv",
    sa + "/quantification_htseq/df_gene_counts_starhtseq.tsv",
]
for s in samples:
    if s != drop:
        files += [sa + "/alignment_htseq/" + s + ".bam",
                  sa + "/alignment_htseq/" + s + ".bam.bai"]
for f in files:
    os.makedirs(os.path.dirname(f), exist_ok=True)
    open(f, "w").write("fake\\n")
"""

FAKE_JAVA = '#!/bin/sh\necho \'openjdk version "17.0.9" 2023-10-17\' >&2\n'
FAKE_DOCKER = "#!/bin/sh\nexit 0\n"
FAKE_SIF = "#!/bin/sh\necho 'apptainer version 1.3.0'\n"


def make_tools(tmp_path: Path, *, java=True, docker=True, singularity=False,
               nextflow=True) -> dict:
    """A restricted PATH holding only the fakes requested plus real git."""
    bindir = tmp_path / "tools"
    bindir.mkdir(exist_ok=True)

    def put(name, body):
        p = bindir / name
        p.write_text(body)
        p.chmod(0o755)

    if nextflow:
        put("nextflow", FAKE_NEXTFLOW)
    if java:
        put("java", FAKE_JAVA)
    if docker:
        put("docker", FAKE_DOCKER)
    if singularity:
        put("apptainer", FAKE_SIF)
    git = shutil.which("git")
    if git and not (bindir / "git").exists():
        (bindir / "git").symlink_to(git)
    return {"PATH": str(bindir)}


GIT = ["git", "-c", "user.name=t", "-c", "user.email=t@t",
       "-c", "commit.gpgsign=false", "-c", "tag.gpgsign=false"]


def make_pipeline(tmp_path: Path, tag: str | None = PIPELINE_TAG,
                  dirty: bool = False) -> Path:
    d = tmp_path / "bj-expression"
    d.mkdir()
    (d / "main.nf").write_text("// fake bj-expression\n")
    (d / "nextflow.config").write_text("manifest { version = '2.0.0' }\n")
    subprocess.run([*GIT, "init", "-q", "-b", "main"], cwd=d, check=True)
    subprocess.run([*GIT, "add", "."], cwd=d, check=True)
    subprocess.run([*GIT, "commit", "-q", "-m", "init"], cwd=d, check=True)
    if tag:
        subprocess.run([*GIT, "tag", tag], cwd=d, check=True)
    if dirty:
        (d / "main.nf").write_text("// edited\n")
    return d


def make_fastqs(tmp_path: Path, names=("s1", "s2")) -> dict:
    d = tmp_path / "fq"
    d.mkdir(exist_ok=True)
    out = {}
    for n in names:
        r1, r2 = d / f"{n}_R1_001.fastq.gz", d / f"{n}_R2_001.fastq.gz"
        r1.write_bytes(b"x")
        r2.write_bytes(b"x")
        out[n] = (r1, r2)
    return out


def write_csv(path: Path, rows, groups=True) -> Path:
    head = "biosampleName,read1,read2" + (",groups" if groups else "")
    lines = [head]
    for name, r1, r2 in rows:
        lines.append(f"{name},{r1},{r2}" + (",g1" if groups else ""))
    path.write_text("\n".join(lines) + "\n")
    return path


def run_cli(args, env=None, expect=0, cwd=None):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, env=env, cwd=cwd,
    )
    if expect is not None:
        assert result.returncode == expect, (
            f"expected exit {expect}, got {result.returncode}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}")
    return result


@pytest.fixture
def setup(tmp_path):
    fq = make_fastqs(tmp_path)
    csv_path = write_csv(tmp_path / "in.csv",
                         [(n, *fq[n]) for n in ("s1", "s2")])
    return {
        "tmp": tmp_path,
        "csv": csv_path,
        "pipeline": make_pipeline(tmp_path),
        "env": make_tools(tmp_path),
        "out": tmp_path / "out",
    }


def base_args(s, *extra):
    return ["--input", s["csv"], "--output", s["out"],
            "--pipeline-dir", s["pipeline"], *extra]


def result_json(out: Path) -> dict:
    return json.loads((out / "result.json").read_text())


# ──────────────────────────────────────────────
# Input CSV validation
# ──────────────────────────────────────────────

def test_valid_csv_is_normalised_to_absolute_paths(tmp_path, monkeypatch):
    mod = _load_module()
    fq = make_fastqs(tmp_path)
    monkeypatch.chdir(tmp_path)
    csv_path = write_csv(tmp_path / "in.csv",
                         [("s1", "fq/s1_R1_001.fastq.gz", "fq/s1_R2_001.fastq.gz")])
    rows = mod.validate_input_csv(csv_path, require_groups=True)
    assert [r["biosampleName"] for r in rows] == ["s1"]
    assert Path(rows[0]["read1"]).is_absolute()
    assert Path(rows[0]["read1"]) == fq["s1"][0].resolve()


def test_missing_required_column_is_named(tmp_path):
    mod = _load_module()
    p = tmp_path / "in.csv"
    p.write_text("biosampleName,read1\ns1,a.fastq.gz\n")
    with pytest.raises(mod.WrapperError, match="read2"):
        mod.validate_input_csv(p, require_groups=False)


def test_header_only_csv_is_rejected(tmp_path):
    mod = _load_module()
    p = tmp_path / "in.csv"
    p.write_text("biosampleName,read1,read2\n")
    with pytest.raises(mod.WrapperError, match="no samples"):
        mod.validate_input_csv(p, require_groups=False)


def test_duplicate_sample_names_are_rejected(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a", "b"))
    p = write_csv(tmp_path / "in.csv",
                  [("dup", *fq["a"]), ("dup", *fq["b"])])
    with pytest.raises(mod.WrapperError, match="dup"):
        mod.validate_input_csv(p, require_groups=True)


@pytest.mark.parametrize("bad", ["has space", "a/b", "..", ".hidden", "-lead",
                                 "semi;colon", ""])
def test_unsafe_sample_names_are_rejected(tmp_path, bad):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv", [(bad, *fq["a"])])
    with pytest.raises(mod.WrapperError):
        mod.validate_input_csv(p, require_groups=True)


def test_missing_local_fastq_is_named(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv",
                  [("a", fq["a"][0], tmp_path / "nope_R2.fastq.gz")])
    with pytest.raises(mod.WrapperError, match="nope_R2"):
        mod.validate_input_csv(p, require_groups=True)


def test_same_file_for_read1_and_read2_is_rejected(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv", [("a", fq["a"][0], fq["a"][0])])
    with pytest.raises(mod.WrapperError, match="same file"):
        mod.validate_input_csv(p, require_groups=True)


def test_s3_paths_are_accepted_and_flagged_remote(tmp_path):
    mod = _load_module()
    p = write_csv(tmp_path / "in.csv",
                  [("a", "s3://bucket/a_R1.fastq.gz", "s3://bucket/a_R2.fastq.gz")])
    rows = mod.validate_input_csv(p, require_groups=True)
    assert rows[0]["read1"] == "s3://bucket/a_R1.fastq.gz"
    assert rows[0]["_remote"] is True


def test_groups_column_is_required_unless_qc_plots_are_skipped(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv", [("a", *fq["a"])], groups=False)
    with pytest.raises(mod.WrapperError, match="groups"):
        mod.validate_input_csv(p, require_groups=True)
    assert mod.validate_input_csv(p, require_groups=False)


# ──────────────────────────────────────────────
# Properties
# ──────────────────────────────────────────────

NAME = st.from_regex(r"[A-Za-z0-9][A-Za-z0-9_-]{0,12}", fullmatch=True)


@settings(deadline=None, max_examples=40)
@given(names=st.lists(NAME, min_size=1, max_size=6, unique=True),
       data=st.data())
@example(names=["a"], data=None)
def test_property_validated_rows_roundtrip_and_ignore_row_order(
        tmp_path_factory, names, data):
    mod = _load_module()
    tmp = tmp_path_factory.mktemp("prop")
    fq = make_fastqs(tmp, tuple(names))
    rows = [(n, *fq[n]) for n in names]
    shuffled = list(rows)
    if data is not None:
        shuffled = data.draw(st.permutations(rows))
    a = mod.validate_input_csv(write_csv(tmp / "a.csv", rows), require_groups=True)
    b = mod.validate_input_csv(write_csv(tmp / "b.csv", shuffled), require_groups=True)
    assert sorted(r["biosampleName"] for r in a) == sorted(names)
    key = lambda r: r["biosampleName"]  # noqa: E731
    assert sorted(a, key=key) == sorted(b, key=key)
    # writing the normalised rows back out and re-reading is stable
    out = tmp / "norm.csv"
    mod.write_normalised_csv(a, out)
    c = mod.validate_input_csv(out, require_groups=True)
    assert sorted(c, key=key) == sorted(a, key=key)


@settings(deadline=None, max_examples=60)
@given(name=st.text(min_size=1, max_size=15))
def test_property_accepted_sample_names_are_path_and_shell_safe(name):
    mod = _load_module()
    if mod.is_safe_sample_name(name):
        assert re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name)
        assert ".." not in name and "/" not in name and " " not in name


# ──────────────────────────────────────────────
# Params and command
# ──────────────────────────────────────────────

def test_params_only_carry_required_values_and_explicit_overrides(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    params = json.loads((setup["out"] / "reproducibility" / "params.json").read_text())
    assert set(params) == {"input_csv", "publish_dir", "timestamp", "genome"}
    assert params["genome"] == "GRCh38"
    for k in ("n_reads", "read_length", "min_reads", "seqtk_sample_seed"):
        assert k not in params, f"{k} would silently override the pipeline default"


def test_explicit_overrides_are_passed_through(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--n-reads", "5000",
                      "--read-length", "75", "--min-reads", "100", "--seed", "7",
                      "--skip-cutadapt", "--genome", "GRCm39",
                      "--max-cpus", "8", "--max-memory", "50.GB"),
            env=setup["env"])
    p = json.loads((setup["out"] / "reproducibility" / "params.json").read_text())
    assert p["n_reads"] == 5000 and p["read_length"] == 75
    assert p["min_reads"] == 100 and p["seqtk_sample_seed"] == 7
    assert p["skip_cutadapt"] is True and p["genome"] == "GRCm39"
    assert p["max_cpus"] == 8 and p["max_memory"] == "50.GB"


@pytest.mark.parametrize("bad", [
    ["--n-reads", "0"], ["--n-reads", "-1"], ["--read-length", "0"],
    ["--min-reads", "-1"], ["--genome", "hg19"], ["--max-cpus", "0"],
    ["--max-memory", "lots"], ["--profile", "podman"],
])
def test_invalid_options_are_rejected_before_anything_runs(setup, bad):
    res = run_cli(base_args(setup, "--accept-bioskryb-terms", *bad),
                  env=setup["env"], expect=None)
    assert res.returncode in (1, 2)
    # the wrapper's / argparse's own message, not a missing-script traceback
    assert "error:" in res.stderr
    assert not (setup["out"] / "result.json").exists()


def test_command_uses_params_file_work_dir_and_pinned_source(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    cmd = (setup["out"] / "reproducibility" / "commands.sh").read_text()
    assert "nextflow run" in cmd and "-params-file" in cmd
    assert str(setup["pipeline"] / "main.nf") in cmd
    assert "-work-dir" in cmd and str(setup["out"]) in cmd
    assert "-resume" not in cmd
    assert "-profile" not in cmd, "docker is the pipeline default"


def test_singularity_profile_and_resume_are_forwarded(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path),
         "env": make_tools(tmp_path, docker=False, singularity=True),
         "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    run_cli(base_args(s, "--accept-bioskryb-terms", "--profile", "singularity",
                      "--resume"), env=s["env"])
    cmd = (s["out"] / "reproducibility" / "commands.sh").read_text()
    assert "-profile singularity" in cmd and "-resume" in cmd


# ──────────────────────────────────────────────
# Pipeline source: never fetched, must be the pinned tag
# ──────────────────────────────────────────────

def test_missing_pipeline_dir_explains_how_to_get_it_and_the_licence(setup):
    res = run_cli(["--input", setup["csv"], "--output", setup["out"],
                   "--pipeline-dir", setup["tmp"] / "nowhere",
                   "--accept-bioskryb-terms"], env=setup["env"], expect=None)
    assert res.returncode != 0
    assert "git clone" in res.stderr and LICENSE_URL in res.stderr


def test_checkout_not_at_the_pinned_tag_is_refused(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path, tag=None),
         "env": make_tools(tmp_path), "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--accept-bioskryb-terms"), env=s["env"], expect=None)
    assert res.returncode != 0
    assert PIPELINE_TAG in res.stderr


def test_dirty_checkout_is_refused_unless_allowed(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path, dirty=True),
         "env": make_tools(tmp_path), "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--accept-bioskryb-terms"), env=s["env"], expect=None)
    assert res.returncode != 0 and "uncommitted" in res.stderr.lower()
    s["out"] = tmp_path / "out2"
    run_cli(base_args(s, "--accept-bioskryb-terms", "--allow-dirty-pipeline"),
            env=s["env"])
    assert result_json(s["out"])["pipeline"]["dirty"] is True


# ──────────────────────────────────────────────
# Licence acknowledgement
# ──────────────────────────────────────────────

def test_running_requires_explicit_acceptance_of_bioskryb_terms(setup):
    log = setup["tmp"] / "nf.log"
    env = {**setup["env"], "FAKE_NF_LOG": str(log)}
    res = run_cli(base_args(setup), env=env, expect=None)
    assert res.returncode != 0
    assert LICENSE_URL in res.stderr
    assert not log.exists(), "nextflow must not be launched without acceptance"


def test_check_mode_does_not_need_acceptance(setup):
    run_cli(base_args(setup, "--check"), env=setup["env"])


# ──────────────────────────────────────────────
# Running (fake nextflow)
# ──────────────────────────────────────────────

def test_successful_run_records_status_outputs_and_samples(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    r = result_json(setup["out"])
    assert r["status"] == "success"
    assert r["terms_accepted"] is True
    assert sorted(r["samples"]) == ["s1", "s2"]
    assert r["samples_missing_from_outputs"] == []
    assert r["outputs"]["multiqc_report"].endswith("multiqc_report.html")
    assert r["outputs"]["salmon_gene_counts"].endswith("df_gene_counts_salmon.tsv")
    assert r["outputs"]["htseq_gene_counts"].endswith("df_gene_counts_starhtseq.tsv")
    assert r["outputs"]["pipeline_metrics_summary"].endswith("pipeline_metrics_summary.csv")
    assert Path(r["outputs"]["results_dir"]).is_dir()
    assert r["pipeline"]["version"] == PIPELINE_TAG


def test_report_has_disclaimer_and_pipeline_licence_notice(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    text = (setup["out"] / "report.md").read_text()
    assert ("ClawBio is a research and educational tool. It is not a medical "
            "device and does not provide clinical diagnoses.") in text
    assert LICENSE_URL in text


def test_samples_dropped_by_the_pipeline_are_called_out(setup):
    env = {**setup["env"], "FAKE_NF_DROP": "s2"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env)
    r = result_json(setup["out"])
    assert r["samples_missing_from_outputs"] == ["s2"]
    text = (setup["out"] / "report.md").read_text().lower()
    assert "s2" in text and "min_reads" in text


def test_nextflow_failure_is_reported_with_stderr_tail(setup):
    env = {**setup["env"], "FAKE_NF_MODE": "fail"}
    res = run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    r = result_json(setup["out"])
    assert r["status"] == "failed" and r["nextflow_exit_code"] == 3
    assert "STARALIGN" in r["stderr_tail"]
    assert "STARALIGN" in res.stderr


def test_exit_zero_without_key_outputs_is_not_success(setup):
    env = {**setup["env"], "FAKE_NF_MODE": "empty"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    r = result_json(setup["out"])
    assert r["status"] == "incomplete"
    assert r["missing_outputs"], "must list what was expected but not produced"


# ──────────────────────────────────────────────
# Preflight and --check
# ──────────────────────────────────────────────

def test_check_validates_without_launching_nextflow(setup):
    log = setup["tmp"] / "nf.log"
    env = {**setup["env"], "FAKE_NF_LOG": str(log)}
    run_cli(base_args(setup, "--check"), env=env)
    r = result_json(setup["out"])
    assert r["status"] == "check_passed"
    assert not (setup["out"] / "results_clawbio").exists()
    # only `nextflow -version` may have been called
    calls = [json.loads(ln) for ln in log.read_text().splitlines()] if log.exists() else []
    assert all(c == ["-version"] for c in calls)


@pytest.mark.parametrize("missing,needle", [
    ({"java": False}, "java"),
    ({"nextflow": False}, "nextflow"),
    ({"docker": False}, "singularity"),
])
def test_missing_tools_are_named_with_a_fix(tmp_path, missing, needle):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path),
         "env": make_tools(tmp_path, **missing), "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--check"), env=s["env"], expect=None)
    assert res.returncode != 0
    assert needle in res.stderr.lower()


def test_singularity_profile_needs_apptainer_or_singularity(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path),
         "env": make_tools(tmp_path, docker=True, singularity=False),
         "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--check", "--profile", "singularity"),
                  env=s["env"], expect=None)
    assert res.returncode != 0
    assert "apptainer" in res.stderr.lower() or "singularity" in res.stderr.lower()


# ──────────────────────────────────────────────
# Output hygiene, reproducibility, demo
# ──────────────────────────────────────────────

def test_existing_results_are_not_overwritten(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    before = (setup["out"] / "result.json").read_bytes()
    res = run_cli(base_args(setup, "--accept-bioskryb-terms"),
                  env=setup["env"], expect=None)
    assert res.returncode != 0
    assert (setup["out"] / "result.json").read_bytes() == before


def test_reproducibility_bundle_is_complete_and_checksums_verify(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    repro = setup["out"] / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256",
                 "params.json", "samplesheet.normalised.csv"):
        assert (repro / name).exists(), name
    for line in (repro / "checksums.sha256").read_text().splitlines():
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((setup["out"] / label).read_bytes()).hexdigest() == digest


def test_writes_only_under_output_dir(setup):
    work = setup["tmp"] / "cwd"
    work.mkdir()
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"], cwd=work)
    assert list(work.iterdir()) == []


def test_demo_needs_no_tools_and_says_nothing_was_run(tmp_path):
    out = tmp_path / "demo"
    run_cli(["--demo", "--output", out], env={"PATH": str(tmp_path)})
    r = result_json(out)
    assert r["status"] == "demo_plan_only"
    assert "plan only" in (out / "report.md").read_text().lower()


def _parse_output_contract(skill_md: Path) -> list[str]:
    text = skill_md.read_text()
    m = re.search(r"## Output Structure\s+```[^\n]*\n(.*?)```", text, re.S)
    if not m:
        return []
    files, stack = [], []
    for line in m.group(1).splitlines()[1:]:
        if not line.strip() or "(optional)" in line:
            continue
        depth = (len(line) - len(line.lstrip("│├└─ "))) // 4 - 1
        name = re.split(r"\s+#", line.lstrip("│├└─ "))[0].strip()
        stack = stack[:depth] + [name.rstrip("/")]
        if not name.endswith("/"):
            files.append("/".join(stack))
    return files


class TestOutputContract:
    def test_documented_outputs_are_produced_by_the_demo(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md has no parseable Output Structure"
        run_cli(["--demo", "--output", tmp_path], env={"PATH": str(tmp_path / "x")})
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, f"promised but not produced: {missing}"
