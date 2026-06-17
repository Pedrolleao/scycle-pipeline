#!/usr/bin/env python3
"""
scycle-pipeline launcher.

Maps MAGs / isolate proteomes to their participation in the sulfur cycle.
Bootstraps the conda environment, detects input type, makes sure the HMM /
BLAST databases are built, writes the discovered samples into
config/config.yaml, then invokes Snakemake.
(Architecture adapted from the sibling Holomicrobiome-ewaste pipeline.)

Inputs are assembled FASTA — .faa proteomes, or .fna assemblies which are
auto-translated with Prodigal. (Raw-FASTQ "Mode B" was removed 2026-06-12 — it
was dead code with no Snakemake rules to run it.)

Usage:
    python run.py [--input PATH] [--mode protein]
                  [--cores N] [--dry-run] [--skip-db-setup]
                  [--prodigal-mode single|meta]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "config" / "config.yaml"
ENV_FILE = ROOT / "envs" / "scycle.yaml"
# Standalone conda env (created on first run from envs/scycle.yaml). Override with
# the SCYCLE_ENV env var to reuse another env with the same deps (e.g. the legacy
# `ewaste-pipeline` env: `SCYCLE_ENV=ewaste-pipeline python run.py …`).
ENV_NAME = os.environ.get("SCYCLE_ENV", "scycle-pipeline")

FASTA_AA_EXT = {".faa", ".fap"}
FASTA_NT_EXT = {".fna", ".fa", ".fasta", ".fas", ".ffn"}
FASTQ_EXT = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

REQUIRED_TOOLS = [
    "snakemake", "hmmscan", "hmmpress", "prodigal", "blastp",
    "makeblastdb", "diamond",
]


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
        sys.exit("error: conda/mamba not found; install miniconda first.")

    # Check whether the env exists.
    out = subprocess.run([conda, "env", "list", "--json"],
                         capture_output=True, text=True, check=True).stdout
    envs = {Path(p).name for p in json.loads(out).get("envs", [])}
    if ENV_NAME not in envs:
        print(f"[run.py] creating conda env '{ENV_NAME}' from {ENV_FILE} "
              "(first run, ~5 min, ~600 MB)…", flush=True)
        subprocess.run([conda, "env", "create", "-f", str(ENV_FILE)], check=True)

    # Re-invoke this script inside the env.
    args = ["conda", "run", "--no-capture-output", "-n", ENV_NAME,
            "python", str(Path(__file__).resolve()), *argv]
    os.execvp(args[0], args)


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
        sys.exit(f"error: {input_path} not found")

    if not paths:
        sys.exit(f"error: no FASTA / FASTQ files in {input_path}")

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
    if len(nt_samples) > 1:
        ans = input(f"Got {len(nt_samples)} nucleotide samples. Treat all as "
                    "[s]ingle isolate genomes or [m]etagenome assemblies? ").strip().lower()
        mode = "meta" if ans.startswith("m") else "single"
        for s in nt_samples:
            s["prodigal_mode"] = mode
        return
    s = nt_samples[0]
    ans = input(f"{s['name']}: [s]ingle isolate or [m]etagenome? ").strip().lower()
    s["prodigal_mode"] = "meta" if ans.startswith("m") else "single"


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
        return
    presence = db_paths_present()
    if not presence["hmm"] or _stale(_HMM_DB):
        why = "missing" if not presence["hmm"] else "stale (targets.yaml newer)"
        print(f"[run.py] building HMM database… ({why})", flush=True)
        subprocess.run([sys.executable, str(ROOT / "workflow" / "scripts" / "build_hmm_db.py"),
                        "--force"], check=True, cwd=ROOT)
    if (not presence["blast_u"] or not presence["blast_g"]
            or _stale(_BLAST_U) or _stale(_BLAST_G)):
        why = ("missing" if not (presence["blast_u"] and presence["blast_g"])
               else "stale (targets.yaml newer)")
        print(f"[run.py] building BLAST databases… ({why})", flush=True)
        subprocess.run([sys.executable, str(ROOT / "workflow" / "scripts" / "build_blast_db.py"),
                        "--force"], check=True, cwd=ROOT)


# ───────────────────────────── config writing ────────────────────────────────

def write_samples_block(samples: list[dict]) -> None:
    """Rewrite the `samples:` block at the bottom of config.yaml.
    Preserves everything else by truncating from the `samples:` marker line."""
    text = CONFIG_PATH.read_text()
    head, _, _ = text.partition("\nsamples:")
    lines = [head, "", "samples:"]
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


# ───────────────────────────── main ──────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="run.py",
        description="scycle-pipeline: map MAGs / proteomes to sulfur-cycle gene participation")
    p.add_argument("--input", type=Path,
                   help="FASTA / FASTQ file or directory")
    p.add_argument("--mode", choices=["auto", "protein"], default="protein",
                   help="Pipeline mode (protein/assembly). Raw-FASTQ 'read' mode was "
                        "removed 2026-06-12; supply .faa proteomes or .fna assemblies.")
    p.add_argument("--cores", type=int, default=os.cpu_count() or 4)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--skip-db-setup", action="store_true")
    p.add_argument("--prodigal-mode", choices=["single", "meta"], default=None,
                   help="Apply to all nt samples (skip the interactive prompt)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Step 1 — env.
    missing = tools_present()
    if missing and not in_env():
        bootstrap_env_and_reexec(sys.argv[1:])  # never returns
    if missing and in_env():
        sys.exit(f"error: missing tools even in env: {missing}")

    # Step 2 — input.
    if args.input is None:
        ans = input("Path to FASTA/FASTQ file or directory: ").strip()
        args.input = Path(ans).expanduser()
    samples = discover_samples(args.input.resolve())

    # Only protein/nucleotide (proteome or assembly) input is supported; raw-reads
    # (FASTQ "Mode B") was removed 2026-06-12 (dead code, no read_mode.smk to run it).
    fastq = [s["name"] for s in samples if s["kind"] == "fastq"]
    if fastq:
        sys.exit(f"error: FASTQ input is not supported (raw-reads mode removed); "
                 f"provide .faa proteomes or .fna assemblies. Offending: {', '.join(fastq)}")

    prompt_prodigal_mode(samples, args.prodigal_mode)

    # Determine pipeline modes touched.
    modes_present = set()
    for s in samples:
        modes_present.add("read" if s["kind"] == "fastq" else "protein")

    # Step 3–5 — databases.
    build_databases(args.skip_db_setup, modes_present)

    # Step 5 — config.
    write_samples_block(samples)
    print(f"[run.py] {len(samples)} sample(s) written to config/config.yaml", flush=True)

    # Step 6 — snakemake.
    cmd = ["snakemake", "--cores", str(args.cores),
           "--configfile", str(CONFIG_PATH),
           "--snakefile", str(ROOT / "workflow" / "Snakefile")]
    if args.dry_run:
        cmd.append("--dry-run")
    print(f"[run.py] $ {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
