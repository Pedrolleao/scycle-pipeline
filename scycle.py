#!/usr/bin/env python3
"""
scycle-pipeline launcher.

Maps MAGs / isolate proteomes to their participation in the sulfur cycle.
Bootstraps the conda environment, detects input type, makes sure the HMM /
BLAST databases are built, writes a per-run config (the base config plus the
discovered samples and the results directory) into the results directory,
then invokes Snakemake with it. The base config is never written.
(Same launcher as the methane and nitrogen sister pipelines.)

Inputs are assembled FASTA — .faa proteomes, or .fna assemblies which are
auto-translated with Prodigal. (Raw-FASTQ "Mode B" was removed 2026-06-12 — it
was dead code with no Snakemake rules to run it.)

Results go to scycle_results/ in this directory (`paths.results_dir` of the config),
or to the directory given with --output.

Terminal output goes through cycle_ui.py (rich, Vee Lab palette). On a terminal the
Snakemake run is shown as one progress bar over its jobs; `--plain`, NO_COLOR or a
non-terminal stdout (logs, CI) give plain text with Snakemake's own output passed
through unchanged. Snakemake's complete log is in .snakemake/log/ either way.

Usage:
    python scycle.py [--input PATH] [--output DIR]
                  [--cores N] [--dry-run] [--skip-db-setup]
                  [--prodigal-mode single|meta] [--plain]

Before Snakemake starts, the launcher refuses what the workflow could not run:
a sample name outside the wildcard pattern of workflow/rules/common.smk, two
inputs with the same stem, a file that is not FASTA, a sequence with characters
outside the FASTA alphabet (diamond blastp hangs on those).
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import cycle_ui as ui

TOOL = "scycle"
SUBTITLE = "Sulfur-cycle gene participation of MAGs and isolate proteomes"

ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "config" / "config.yaml"
ENV_FILE = ROOT / "envs" / "scycle.yaml"
# Conda env, shared with the nitrogen (ncycle-pipeline) sister pipeline: both tools
# have the same dependencies, so one env named `cycle-pipeline` serves both. It is
# created on first run from envs/scycle.yaml if it does not exist yet. Override
# with the SCYCLE_ENV env var to use a different env with the same deps.
ENV_NAME = os.environ.get("SCYCLE_ENV", "cycle-pipeline")

FASTA_AA_EXT = {".faa", ".fap"}
FASTA_NT_EXT = {".fna", ".fa", ".fasta", ".fas", ".ffn"}
FASTQ_EXT = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

# A sample name is the {sample} wildcard of every rule. Same pattern as the
# `wildcard_constraints` of workflow/rules/common.smk — change both together.
SAMPLE_NAME_RE = re.compile(r"[A-Za-z0-9_.\-]+")
SAMPLE_NAME_RULE = "letters, digits, '_', '.', '-' only (no spaces)"
# Residues a FASTA sequence line may hold (amino acids or nucleotides, IUPAC codes,
# '*' stop, '-' / '.' gaps). Anything else makes diamond blastp hang.
FASTA_ALPHABET = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*-.")

REQUIRED_TOOLS = [
    "snakemake", "hmmscan", "hmmpress", "prodigal", "blastp",
    "makeblastdb", "diamond",
]

N_STAGES = 4


# ───────────────────────────── environment bootstrap ─────────────────────────

def in_env() -> bool:
    return os.environ.get("CONDA_DEFAULT_ENV") == ENV_NAME


def tools_present() -> list[str]:
    return [t for t in REQUIRED_TOOLS if shutil.which(t) is None]


def bootstrap_env_and_reexec(argv: list[str]) -> None:
    """If we are not inside the target conda env (ENV_NAME), create it (if
    needed) and re-invoke this script via `conda run -n …`. Never returns."""
    conda = shutil.which("conda") or shutil.which("mamba")
    if conda is None:
        ui.err("conda/mamba not found; install miniconda first.")
        sys.exit(1)

    # Check whether the env exists.
    out = subprocess.run([conda, "env", "list", "--json"],
                         capture_output=True, text=True, check=True).stdout
    envs = {Path(p).name for p in json.loads(out).get("envs", [])}
    if ENV_NAME not in envs:
        ui.info(f"creating conda env '{ENV_NAME}' from {ENV_FILE} (first run, ~5 min, ~600 MB)…")
        subprocess.run([conda, "env", "create", "-f", str(ENV_FILE)], check=True)
    else:
        ui.info(f"re-running inside conda env '{ENV_NAME}'")

    # Re-invoke this script inside the env.
    args = [conda, "run", "--no-capture-output", "-n", ENV_NAME,
            "python", str(Path(__file__).resolve()), *argv]
    os.execvp(args[0], args)


def tool_version() -> str:
    """Short git revision of the checkout, '' when not a git checkout."""
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "describe", "--tags", "--always", "--dirty"],
                             capture_output=True, text=True, timeout=5)
        rev = out.stdout.strip()
        return f"rev {rev}" if out.returncode == 0 and rev else ""
    except (OSError, subprocess.SubprocessError):
        return ""


# ───────────────────────────── input detection ───────────────────────────────

def looks_like_fastq(path: Path) -> bool:
    suffixes = "".join(path.suffixes[-2:]).lower()
    return any(suffixes.endswith(e) for e in FASTQ_EXT)


def sample_fasta(path: Path, sample_size: int = 3000) -> tuple[bool, str]:
    """(first non-blank line is a '>' header, up to `sample_size` residues of the
    first records). One read serves the type detection and the sanity check."""
    buf: list[str] = []
    has_header = False
    first = True
    with open(path, errors="replace") as fh:
        for line in fh:
            if first and line.strip():
                has_header = line.startswith(">")
                first = False
            if line.startswith(">"):
                continue
            buf.append(line.strip())
            if sum(len(s) for s in buf) >= sample_size:
                break
    return has_header, "".join(buf)[:sample_size]


def fasta_is_nucleotide(path: Path, sample_size: int = 3000) -> bool:
    """Read up to `sample_size` non-header residues; nucleotide if >85% are
    in {A,C,G,T,N,U} (case-insensitive)."""
    _, seq = sample_fasta(path, sample_size)
    seq = seq.upper()
    if not seq:
        return False
    nt = sum(1 for c in seq if c in "ACGTUN")
    return nt / len(seq) > 0.85


def fasta_problem(path: Path) -> str | None:
    """Why the workflow could not run this file, or None. Looks at the first
    records only (what sample_fasta reads): no '>' header, no sequence, or a
    character outside the FASTA alphabet."""
    has_header, seq = sample_fasta(path)
    if not has_header:
        return "not a FASTA file (the first line is not a '>' header)"
    if not seq:
        return "no sequence after the first header"
    bad = next((c for c in seq if c not in FASTA_ALPHABET), None)
    if bad is not None:
        return f"illegal character {bad!r} in the sequence (diamond blastp would hang on it)"
    return None


def validate_samples(samples: list[dict], check_sequences: bool = True) -> None:
    """Refuse what the workflow cannot run; every problem is listed, then exit 1.
    Names are not fixed silently: the sample name is the row label of every output."""
    problems: list[str] = []
    by_name: dict[str, list[str]] = {}
    for s in samples:
        by_name.setdefault(s["name"], []).append(s["path"])
        if not SAMPLE_NAME_RE.fullmatch(s["name"]):
            problems.append(f"{s['path']}: sample name '{s['name']}' is not allowed — "
                            f"{SAMPLE_NAME_RULE}; rename the file")
    for name, paths in by_name.items():
        if len(paths) > 1:
            problems.append(f"two inputs would be one sample '{name}': "
                            + " and ".join(paths) + "; rename or remove one")
    if check_sequences:
        for s in samples:
            if s["kind"] == "fastq":
                continue
            why = fasta_problem(Path(s["path"]))
            if why:
                problems.append(f"{s['path']}: {why}")
    if problems:
        for msg in problems:
            ui.err(msg)
        ui.err(f"{len(problems)} input problem(s); nothing was run")
        sys.exit(1)


def discover_samples(input_path: Path) -> list[dict]:
    """Scan a file or directory and return a list of sample dicts:
    { name, path, kind } where kind ∈ {"protein", "nucleotide", "fastq"}."""
    paths: list[Path] = []
    if input_path.is_file():
        paths = [input_path]
    elif input_path.is_dir():
        for p in sorted(input_path.iterdir()):
            if not p.is_file():
                continue
            suffixes = "".join(p.suffixes[-2:]).lower()
            if (any(suffixes.endswith(e) for e in FASTQ_EXT) or
                p.suffix.lower() in FASTA_AA_EXT | FASTA_NT_EXT):
                paths.append(p)
    else:
        ui.err(f"{input_path} not found")
        sys.exit(1)

    if not paths:
        ui.err(f"no FASTA / FASTQ files in {input_path}")
        sys.exit(1)

    out = []
    seen: dict[str, list[Path]] = {}
    for p in paths:
        if looks_like_fastq(p):
            # group paired-end by base name (strip _R1/_R2)
            base = p.name
            for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
                if base.lower().endswith(ext):
                    base = base[: -len(ext)]
                    break
            for tag in ("_R1", "_R2", "_1", "_2"):
                if base.endswith(tag):
                    base = base[: -len(tag)]
                    break
            seen.setdefault(base, []).append(p)
            continue

        kind = "nucleotide" if fasta_is_nucleotide(p) else "protein"
        name = p.stem
        out.append({"name": name, "path": str(p.resolve()), "kind": kind})

    for base, files in seen.items():
        files = sorted(files)
        r1 = next((f for f in files if any(t in f.name for t in ("_R1", "_1."))),
                  files[0])
        r2 = next((f for f in files if any(t in f.name for t in ("_R2", "_2."))),
                  None)
        sample = {"name": base, "kind": "fastq",
                  "path": str(r1.resolve())}
        if r2 is not None and r2 != r1:
            sample["path2"] = str(r2.resolve())
        out.append(sample)
    return out


PRODIGAL_CHOICES = {"single": "isolate genome", "meta": "metagenome assembly"}


def prompt_prodigal_mode(samples: list[dict],
                         override: str | None) -> None:
    """For nucleotide samples, fill in `prodigal_mode` (single|meta).
    `override` from CLI applies to all; otherwise prompt the user."""
    nt_samples = [s for s in samples if s["kind"] == "nucleotide"]
    if not nt_samples:
        return
    if override:
        for s in nt_samples:
            s["prodigal_mode"] = override
        return
    if not ui.can_ask():
        # Batch job (cluster, CI): nothing to answer on. Use the default and say so.
        for s in nt_samples:
            s["prodigal_mode"] = "single"
        ui.warn(f"{len(nt_samples)} nucleotide sample(s) and no --prodigal-mode: no terminal to "
                "ask on, so Prodigal runs in 'single' mode (isolate genomes). For metagenome "
                "assemblies pass --prodigal-mode meta.")
        return
    if len(nt_samples) > 1:
        mode = ui.choose(f"Got {len(nt_samples)} nucleotide samples. Treat all as "
                         "single isolate genomes or metagenome assemblies?",
                         PRODIGAL_CHOICES, default="single")
        for s in nt_samples:
            s["prodigal_mode"] = mode
        return
    s = nt_samples[0]
    s["prodigal_mode"] = ui.choose(f"{s['name']}: single isolate or metagenome?",
                                   PRODIGAL_CHOICES, default="single")


def samples_table(samples: list[dict]) -> None:
    rows = []
    for s in samples:
        p = Path(s["path"])
        try:
            size = ui.fmt_size(p.stat().st_size)
        except OSError:
            size = "?"
        prodigal = s.get("prodigal_mode", "—" if s["kind"] == "protein" else "?")
        rows.append([s["name"], s["kind"], size, prodigal])
    ui.table(f"{len(samples)} sample(s)", ["sample", "input", "size", "Prodigal"], rows,
             justify={"input": "left", "Prodigal": "left"})


# ───────────────────────────── DB freshness ──────────────────────────────────

_HMM_DB   = ROOT / "resources" / "hmm" / "scycle_targets.hmm.h3i"
_BLAST_U  = ROOT / "resources" / "blast_db" / "unstable_refs.phr"
_BLAST_G  = ROOT / "resources" / "blast_db" / "blast_gated_refs.phr"
_TARGETS  = ROOT / "config" / "targets.yaml"


def db_paths_present() -> dict[str, bool]:
    return {
        "hmm":      _HMM_DB.exists(),
        "blast_u":  _BLAST_U.exists(),
        "blast_g":  _BLAST_G.exists(),
    }


def _stale(db: Path) -> bool:
    """A built DB is stale if config/targets.yaml (the marker/seed/threshold
    source of truth) was edited after the DB was last built. The DB build runs
    OUTSIDE the Snakemake DAG, so without this check an edit to targets.yaml
    (a new KO, a seed change, a TC override) would silently leave the old DB in
    place. The custom-HMM splice and the seed cache are downstream of
    targets.yaml, so its mtime is the conservative trigger."""
    return db.exists() and _TARGETS.exists() and _TARGETS.stat().st_mtime > db.stat().st_mtime


def build_databases(skip: bool) -> None:
    if skip:
        ui.info("database check skipped (--skip-db-setup)")
        return
    presence = db_paths_present()
    if not presence["hmm"] or _stale(_HMM_DB):
        why = "missing" if not presence["hmm"] else "stale (targets.yaml newer)"
        ui.info(f"building the HMM database ({why}) — the build script's output follows")
        subprocess.run([sys.executable, str(ROOT / "workflow" / "scripts" / "build_hmm_db.py"),
                        "--force"], check=True, cwd=ROOT)
        ui.ok(f"HMM database built: {_HMM_DB.parent}")
    else:
        ui.ok(f"HMM database up to date: {_HMM_DB.parent}")
    if (not presence["blast_u"] or not presence["blast_g"]
            or _stale(_BLAST_U) or _stale(_BLAST_G)):
        why = ("missing" if not (presence["blast_u"] and presence["blast_g"])
               else "stale (targets.yaml newer)")
        ui.info(f"building the BLAST databases ({why}) — the build script's output follows")
        subprocess.run([sys.executable, str(ROOT / "workflow" / "scripts" / "build_blast_db.py"),
                        "--force"], check=True, cwd=ROOT)
        ui.ok(f"BLAST databases built: {_BLAST_U.parent}")
    else:
        ui.ok(f"BLAST databases up to date: {_BLAST_U.parent}")


# ───────────────────────────── config writing ────────────────────────────────

RUN_CONFIG_NAME = f"{TOOL}_run_config.yaml"


def write_run_config(samples: list[dict], results_dir: Path, version: str,
                     explicit_output: bool, dry_run: bool = False) -> Path:
    """Write <results_dir>/scycle_run_config.yaml: the base config (everything above
    its `samples:` line) with the discovered `samples:` block and — when --output was
    given — `results_dir:`; the tool version, date and command line are comments at
    the top. Snakemake is given this file; the base config is only read, so a tracked
    config stays clean and two runs with different results directories do not clobber
    each other. The file doubles as the provenance record of the run, so a dry run
    writes a temporary file instead of overwriting the record of the last real run."""
    text = CONFIG_PATH.read_text()
    head, _, _ = text.partition("\nsamples:")
    lines = [
        f"# {RUN_CONFIG_NAME} — written by {TOOL}.py, do not edit: the next run overwrites it.",
        f"# base config: {CONFIG_PATH}",
        f"# results: {results_dir}",
        f"# tool: {TOOL} {version or '(no git revision)'}",
        f"# date: {datetime.datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"# command: {' '.join(sys.argv)}",
        "",
        head.rstrip("\n"),
        "",
    ]
    if explicit_output:
        # Without --output, paths.results_dir of the base config decides (relative to
        # the pipeline directory, as before), so existing results are not re-keyed.
        lines += ["# Results directory of this run (--output; overrides paths.results_dir).",
                  f"results_dir: {results_dir}", ""]
    lines.append("samples:" if samples else "samples: {}")
    for s in samples:
        lines.append(f"  {s['name']}:")
        for k, v in s.items():
            if k != "name":
                lines.append(f"    {k}: {v}")
    if dry_run:
        fd, name = tempfile.mkstemp(prefix=f"{TOOL}_run_config.", suffix=".yaml")
        os.close(fd)
        out = Path(name)
    else:
        results_dir.mkdir(parents=True, exist_ok=True)
        out = results_dir / RUN_CONFIG_NAME
    out.write_text("\n".join(lines) + "\n")
    return out


# ───────────────────────────── snakemake ─────────────────────────────────────

_RE_STEPS = re.compile(r"^(\d+) of (\d+) steps \((\d+)%\) done")
_RE_RULE = re.compile(r"^(?:local)?rule (\w+):")
_RE_FINISHED = re.compile(r"^Finished jobid: \d+ \(Rule: (\w+)\)")
_RE_STAT_ROW = re.compile(r"^(\w+)\s+(\d+)$")
# Snakemake 9 prints "Complete log(s):", older versions "Complete log:".
_RE_LOG = re.compile(r"Complete log(?:\(s\)| file)?:\s*(\S+)")
_RE_ERROR = re.compile(r"^(Error in rule \w+|\w*Error\b|\w*Exception\b|Shutting down|Exiting because)")
_RE_ERROR_IN_RULE = re.compile(r"^Error in rule \w+")
_RE_JOB_LOG = re.compile(r"^\s*log:\s*(\S+)")
# Snakemake's job header lines: no cause in them, dropped from the context shown.
_RE_CHATTER = re.compile(r"^(\s+(input|output|log|jobid|reason|wildcards|threads|resources):|"
                         r"(local)?rule \w+:|\[\w{3} \w{3}|Select jobs|Execute \d|\s*$)")


def failure_report(error_block: list[str], context: list[str]) -> tuple[str, str]:
    """(head line, panel body) for a failed run. The cause of a failed shell job is in
    the job's log file, not in Snakemake's block (which repeats the command), so the
    body starts with the last lines of that log when the rule has one; for a rule
    without `log:` the tool's own output arrived on the pipe just before the block and
    the last lines of that are shown instead. Snakemake's description of the job
    follows, up to the `shell:` line, once (Snakemake 9 prints the block twice)."""
    if not error_block:
        return "", ""
    head = next((l for l in error_block if _RE_ERROR_IN_RULE.match(l)), error_block[0])
    desc: list[str] = []
    for line in error_block[error_block.index(head):]:
        if line.strip().startswith("shell:") or line.startswith(("Shutting down", "Exiting because")):
            break
        if desc and _RE_ERROR_IN_RULE.match(line):
            break
        if line.strip() != "message: None":
            desc.append(line.rstrip())
    if any("(command exited with non-zero exit code)" in l for l in error_block):
        desc.append("    (command exited with non-zero exit code)")

    body: list[str] = []
    m = next((_RE_JOB_LOG.match(l) for l in error_block if _RE_JOB_LOG.match(l)), None)
    job_log = Path(m.group(1)) if m else None
    if job_log is not None and job_log.is_file() and job_log.stat().st_size > 0:
        tail = job_log.read_text(errors="replace").splitlines()[-15:]
        body += [f"job log — {job_log}:", *tail, ""]
    else:
        ctx = [l.rstrip() for l in context if not _RE_CHATTER.match(l)]
        if ctx:
            body += ["last output before the error:", *ctx[-25:], ""]
        if job_log is not None:
            body.append(f"(the job log {job_log} is empty)")
    body += ["Snakemake:", *desc]
    return head, "\n".join(body)


def run_snakemake_live(cmd: list[str]) -> int:
    """Run Snakemake with its output captured, showing one progress bar over the jobs.
    The lines are parsed for the job plan ("Job stats"), rule starts/finishes and the
    "N of M steps" counter; everything else stays in Snakemake's own log. On failure the
    "Error in rule" block is shown. Returns Snakemake's exit code."""
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    proc = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1, env=env)
    lines: list[str] = []
    stats: dict[str, int] = {}
    in_stats = False
    plan_shown = False
    running: list[str] = []
    error_block: list[str] = []
    context: list[str] = []      # pipe lines just before the error block
    in_error = False
    nothing_to_do = False
    log_path = ""
    done = 0

    def description() -> str:
        if not running:
            return "waiting for the next job"
        names = sorted(set(running))
        shown = ", ".join(names[:3]) + (f" +{len(names) - 3}" if len(names) > 3 else "")
        return f"running {shown}"

    try:
        with ui.progress() as prog:
            task = prog.add_task("building the DAG of jobs", total=None)
            assert proc.stdout is not None
            for raw in proc.stdout:
                line = raw.rstrip("\n")
                lines.append(line)

                m = _RE_LOG.search(line)
                if m:
                    log_path = m.group(1)
                    continue
                if in_error and _RE_STEPS.match(line):
                    # Snakemake went on: the "error" was a tool's own message, not fatal.
                    in_error, error_block, context = False, [], []
                elif in_error:
                    if not line.startswith(("Shutting down", "Exiting because")):
                        error_block.append(line)
                    continue
                if _RE_ERROR.match(line):
                    in_error = True
                    context = lines[-26:-1]
                    error_block.append(line)
                    continue

                if line.startswith("Job stats:"):
                    in_stats = True
                    continue
                if in_stats:
                    m = _RE_STAT_ROW.match(line)
                    if m and m.group(1) not in ("job",):
                        stats[m.group(1)] = int(m.group(2))
                    if m and m.group(1) == "total":
                        in_stats = False
                        n_jobs = stats.pop("total", 0) - (1 if "all" in stats else 0)
                        stats.pop("all", None)
                        prog.update(task, total=max(n_jobs, 1), description=description())
                        if not plan_shown:
                            rows = [[r, c] for r, c in sorted(stats.items(), key=lambda kv: (-kv[1], kv[0]))]
                            prog.console.print()
                            ui.table(f"Plan: {n_jobs} job(s)", ["rule", "jobs"], rows)
                            plan_shown = True
                    continue

                m = _RE_RULE.match(line)
                if m:
                    if m.group(1) != "all":
                        running.append(m.group(1))
                    prog.update(task, description=description())
                    continue
                m = _RE_FINISHED.match(line)
                if m:
                    if m.group(1) in running:
                        running.remove(m.group(1))
                    prog.update(task, description=description())
                    continue
                m = _RE_STEPS.match(line)
                if m:
                    done, total = int(m.group(1)), int(m.group(2))
                    # `all` is the last step and produces nothing; count the real jobs.
                    real_total = total - 1 if total > 1 else total
                    prog.update(task, completed=min(done, real_total), total=real_total,
                                description=description())
                    continue
                if line.startswith("Nothing to be done"):
                    nothing_to_do = True
                    prog.update(task, total=1, completed=1, description="nothing to do")
                    continue
            rc = proc.wait()
            if rc == 0 and not nothing_to_do:
                prog.update(task, completed=prog.tasks[0].total, description="all jobs finished")
    except KeyboardInterrupt:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        ui.warn("interrupted — Snakemake was stopped; re-run to resume from the finished jobs")
        return 130

    if nothing_to_do:
        ui.ok("nothing to do — the results are already up to date")
    if rc != 0:
        if error_block:
            head, tail = failure_report(error_block, context)
        else:
            head, tail = f"Snakemake exited with status {rc}", "\n".join(lines[-25:])
        ui.error_panel("Pipeline failed", head, tail,
                       hint=f"complete log: {log_path or '.snakemake/log/'}")
    return rc


def run_snakemake(cmd: list[str], dry_run: bool) -> int:
    """Live rendering on a terminal; pass-through of Snakemake's own output for
    --dry-run, --plain, NO_COLOR or when stdout is not a terminal."""
    if dry_run or not ui.is_live():
        ui.info("Snakemake's output follows" + (" (dry run)" if dry_run else ""))
        return subprocess.run(cmd, cwd=ROOT).returncode
    return run_snakemake_live(cmd)


def results_summary(results_dir: Path, samples: list[dict]) -> None:
    files = []
    for label, rel in (("report", "scycle_report.html"),
                       ("matrix", "scycle_matrix.tsv"),
                       ("heatmap", "scycle_heatmap.svg"),
                       ("figures", "figures/")):
        if (results_dir / rel.rstrip("/")).exists():
            files.append((label, str(results_dir / rel)))
    files.append(("per sample", f"{results_dir}/<sample>/report/"))
    ui.results_panel("Results", f"{len(samples)} sample(s) mapped to the sulfur cycle", files,
                     hint="open the report first; gap_analysis.txt in each sample's report/ explains the calls")


# ───────────────────────────── main ──────────────────────────────────────────

def configured_results_dir() -> Path:
    """`paths.results_dir` of the config file, as an absolute path."""
    m = re.search(r"^\s*results_dir:\s*(\S+)", CONFIG_PATH.read_text(), re.M)
    return (ROOT / (m.group(1).strip("'\"") if m else "scycle_results")).resolve()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="scycle.py",
        description="scycle-pipeline: map MAGs / proteomes to sulfur-cycle gene participation")
    p.add_argument("--input", type=Path,
                   help="FASTA / FASTQ file or directory")
    p.add_argument("--output", type=Path, default=None, metavar="DIR",
                   help="Directory for the results (default: scycle_results/ next to "
                        "this script, as set by paths.results_dir in the config)")
    # Only the protein/assembly mode exists (raw-FASTQ 'read' mode was removed
    # 2026-06-12); the flag is kept for old command lines and hidden from --help.
    p.add_argument("--mode", choices=["auto", "protein"], default="protein",
                   help=argparse.SUPPRESS)
    # Skips the FASTA sanity check (names and duplicates are still checked): lets a
    # test feed a bad file to the workflow to see how a failing job is reported.
    p.add_argument("--skip-validation", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--cores", type=int, default=os.cpu_count() or 4)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--skip-db-setup", action="store_true")
    p.add_argument("--prodigal-mode", choices=["single", "meta"], default=None,
                   help="Apply to all nt samples (skip the interactive prompt)")
    p.add_argument("--plain", action="store_true",
                   help="Plain text output: no colours, no progress bar, Snakemake's own "
                        "output passed through (also the default when stdout is not a terminal)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ui.init(TOOL, plain=args.plain)

    # Step 1 — env.
    missing = tools_present()
    if missing and not in_env():
        bootstrap_env_and_reexec(sys.argv[1:])  # never returns
    if missing and in_env():
        ui.err(f"missing tools even in env '{ENV_NAME}': {', '.join(missing)}")
        sys.exit(1)

    version = tool_version()
    ui.banner("run" if not args.dry_run else "dry run", version, SUBTITLE)

    # --output is relative to where the user is; Snakemake runs from ROOT.
    results_dir = (args.output.expanduser().resolve() if args.output is not None
                   else configured_results_dir())

    with ui.stage(1, N_STAGES, "Environment"):
        ui.ok(f"conda env '{os.environ.get('CONDA_DEFAULT_ENV', ENV_NAME)}' — "
              f"{len(REQUIRED_TOOLS)} required tools found")
        ui.kv("config", CONFIG_PATH)
        ui.kv("results", results_dir)
        ui.kv("cores", args.cores)

    # Step 2 — input.
    with ui.stage(2, N_STAGES, "Input"):
        if args.input is None:
            if not ui.can_ask():
                ui.err("no --input given and no terminal to ask on (batch job?): "
                       "run  python scycle.py --input PATH [--prodigal-mode single|meta]")
                sys.exit(2)
            try:
                ans = ui.ask("Path to a FASTA file or directory:")
            except (EOFError, KeyboardInterrupt):
                print()
                sys.exit(130)
            args.input = Path(ans).expanduser()
        ui.kv("input", args.input.resolve())
        samples = discover_samples(args.input.resolve())

        # Only protein/nucleotide (proteome or assembly) input is supported; raw-reads
        # (FASTQ "Mode B") was removed 2026-06-12 (dead code, no read_mode.smk to run it).
        fastq = [s["name"] for s in samples if s["kind"] == "fastq"]
        if fastq:
            ui.err("FASTQ input is not supported (raw-reads mode removed); "
                   f"provide .faa proteomes or .fna assemblies. Offending: {', '.join(fastq)}")
            sys.exit(1)

        validate_samples(samples, check_sequences=not args.skip_validation)
        prompt_prodigal_mode(samples, args.prodigal_mode)
        samples_table(samples)

    # Step 3 — databases.
    with ui.stage(3, N_STAGES, "Databases"):
        build_databases(args.skip_db_setup)

    # Step 4 — the per-run config, then Snakemake.
    run_config = write_run_config(samples, results_dir, version, args.output is not None,
                                  dry_run=args.dry_run)
    cmd = ["snakemake", "--cores", str(args.cores),
           "--configfile", str(run_config),
           "--snakefile", str(ROOT / "workflow" / "Snakefile")]
    if args.dry_run:
        cmd.append("--dry-run")

    with ui.stage(4, N_STAGES, "Pipeline"):
        ui.ok(f"{len(samples)} sample(s) written to {run_config}")
        ui.kv("command", " ".join(cmd))
        rc = run_snakemake(cmd, args.dry_run)
    if rc != 0:
        sys.exit(rc)
    if not args.dry_run:
        results_summary(results_dir, samples)


if __name__ == "__main__":
    main()
