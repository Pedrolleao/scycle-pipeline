#!/usr/bin/env python3
"""
scycle-pipeline launcher.

Maps MAGs / isolate proteomes to their participation in the sulfur cycle.
Bootstraps the conda environment, detects input type, makes sure the HMM /
BLAST databases are built, writes the discovered samples into
config/config.yaml, then invokes Snakemake.
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
    python scycle.py [--input PATH] [--output DIR] [--mode protein]
                  [--cores N] [--dry-run] [--skip-db-setup]
                  [--prodigal-mode single|meta] [--plain]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
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
    args = ["conda", "run", "--no-capture-output", "-n", ENV_NAME,
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


def fasta_is_nucleotide(path: Path, sample_size: int = 3000) -> bool:
    """Read up to `sample_size` non-header residues; nucleotide if >85% are
    in {A,C,G,T,N,U} (case-insensitive)."""
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            buf.append(line.strip())
            if sum(len(s) for s in buf) >= sample_size:
                break
    seq = "".join(buf)[:sample_size].upper()
    if not seq:
        return False
    nt = sum(1 for c in seq if c in "ACGTUN")
    return nt / len(seq) > 0.85


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


def build_databases(skip: bool, modes_present: set[str]) -> None:
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

def write_samples_block(samples: list[dict]) -> None:
    """Rewrite the `samples:` block at the bottom of config.yaml.
    Preserves everything else by truncating from the `samples:` marker line."""
    text = CONFIG_PATH.read_text()
    head, _, _ = text.partition("\nsamples:")
    lines = [head.rstrip("\n"), "", "samples:"]
    if not samples:
        lines[-1] = "samples: {}"
    else:
        for s in samples:
            lines.append(f"  {s['name']}:")
            for k, v in s.items():
                if k == "name":
                    continue
                lines.append(f"    {k}: {v}")
    CONFIG_PATH.write_text("\n".join(lines) + "\n")


# ───────────────────────────── snakemake ─────────────────────────────────────

_RE_STEPS = re.compile(r"^(\d+) of (\d+) steps \((\d+)%\) done")
_RE_RULE = re.compile(r"^(?:local)?rule (\w+):")
_RE_FINISHED = re.compile(r"^Finished jobid: \d+ \(Rule: (\w+)\)")
_RE_STAT_ROW = re.compile(r"^(\w+)\s+(\d+)$")
_RE_LOG = re.compile(r"Complete log(?: file)?:\s*(\S+)")
_RE_ERROR = re.compile(r"^(Error in rule \w+|\w*Error\b|\w*Exception\b|Shutting down|Exiting because)")


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
                    in_error, error_block = False, []
                elif in_error:
                    if not line.startswith(("Shutting down", "Exiting because")):
                        error_block.append(line)
                    continue
                if _RE_ERROR.match(line):
                    in_error = True
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
        head = error_block[0] if error_block else f"Snakemake exited with status {rc}"
        tail = "\n".join(error_block[1:] if error_block else lines[-25:])
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
    p.add_argument("--mode", choices=["auto", "protein"], default="protein",
                   help="Pipeline mode (protein/assembly). Raw-FASTQ 'read' mode was "
                        "removed 2026-06-12; supply .faa proteomes or .fna assemblies.")
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

    ui.banner("run" if not args.dry_run else "dry run", tool_version(), SUBTITLE)

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

        prompt_prodigal_mode(samples, args.prodigal_mode)
        samples_table(samples)

    # Determine pipeline modes touched.
    modes_present = set()
    for s in samples:
        modes_present.add("read" if s["kind"] == "fastq" else "protein")

    # Step 3–5 — databases.
    with ui.stage(3, N_STAGES, "Databases"):
        build_databases(args.skip_db_setup, modes_present)

    # Step 5 — config.
    write_samples_block(samples)

    # Step 6 — snakemake.
    cmd = ["snakemake", "--cores", str(args.cores),
           "--configfile", str(CONFIG_PATH),
           "--snakefile", str(ROOT / "workflow" / "Snakefile")]
    if args.output is not None:
        cmd += ["--config", f"results_dir={results_dir}"]
    if args.dry_run:
        cmd.append("--dry-run")

    with ui.stage(4, N_STAGES, "Pipeline"):
        ui.ok(f"{len(samples)} sample(s) written to {CONFIG_PATH.relative_to(ROOT) if CONFIG_PATH.is_relative_to(ROOT) else CONFIG_PATH}")
        ui.kv("command", " ".join(cmd))
        rc = run_snakemake(cmd, args.dry_run)
    if rc != 0:
        sys.exit(rc)
    if not args.dry_run:
        results_summary(results_dir, samples)


if __name__ == "__main__":
    main()
