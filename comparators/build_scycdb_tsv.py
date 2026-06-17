#!/usr/bin/env python3
"""
build_scycdb_tsv.py — run SCycDB (Yu et al. 2020) on the SP1 panel and export the
normalized TSV the benchmark harness consumes (columns: sample, family, count).

SCycDB is the sulfur-cycle domain database analogous to NCycDB for nitrogen
(same authors, qichao1984; same file format). Faithful to SCycDB's own profiler
(SCycDB_FunctionProfiler.PL), protein mode, DEFAULT settings:
    diamond makedb --in data/SCycDB_2020Mar.faa --db data/SCycDB_2020Mar   (once)
    diamond blastp -k 1 -e 1e-4 -d data/SCycDB_2020Mar -q <proteome> -o <hits>
(nucleotide inputs use `diamond blastx`, matching the profiler's seqtype handling).
Best-hit subject id -> gene family via data/id2gene.2020Mar.map; a family's `count`
for a genome = number of query proteins whose top hit lands in that family
(present iff > 0). Mirrors ../../comparators? no — comparators/build_ncycdb_tsv.py
in the nitrogen tool.

SCycDB provenance: github.com/qichao1984/SCycDB, SCycDB_2020Mar.faa (911,805 seqs;
579,056 family-mapped + homolog decoys), id2gene.2020Mar.map (207 gene families).
DB + per-genome hits cached under comparators/SCyc/ + comparators/scycdb_out/.
"""
from __future__ import annotations
import csv, subprocess, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCYC = HERE / "SCyc" / "data"
DB = SCYC / "SCycDB_2020Mar"                  # diamond db basename (built from .faa)
FAA = SCYC / "SCycDB_2020Mar.faa"
ID2GENE = SCYC / "id2gene.2020Mar.map"
ROOT = HERE.parent                            # scycle-pipeline/
PANEL = ROOT.parent / "sp1_panel"             # Sulfur_Cycle/sp1_panel/*.faa (+ 1 .fna)
OUTDIR = HERE / "scycdb_out"
OUT_TSV = ROOT / "validation" / "benchmark" / "scycdb.tsv"
EVALUE = "0.0001"                             # SCycDB_FunctionProfiler.PL default (-e 1e-4)
THREADS = "8"


def ensure_db() -> None:
    if not DB.with_suffix(".dmnd").exists():
        print(f"[scycdb] building DIAMOND DB from {FAA.name} (~one-time)…", file=sys.stderr)
        subprocess.run(["diamond", "makedb", "--in", str(FAA), "--db", str(DB), "--quiet"],
                       check=True)


def load_id2gene() -> dict[str, str]:
    m = {}
    with open(ID2GENE) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                m[p[0]] = p[1]
    return m


def panel_inputs(proteome_dir: Path) -> list[tuple[str, Path]]:
    """(genome_stem, path) for every proteome (.faa) or assembly (.fna) in the dir."""
    out = []
    for p in sorted(proteome_dir.iterdir()):
        if p.suffix.lower() in (".faa", ".fna"):
            out.append((p.stem, p))
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Run SCycDB on a proteome panel -> normalized TSV.")
    ap.add_argument("--proteome-dir", type=Path, default=PANEL,
                    help="dir of <genome>.faa/.fna (default: SP1 panel)")
    ap.add_argument("--out", type=Path, default=OUT_TSV, help="output normalized TSV")
    ap.add_argument("--outdir", type=Path, default=OUTDIR, help="per-genome DIAMOND-hits cache dir")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    ensure_db()
    id2gene = load_id2gene()
    inputs = panel_inputs(args.proteome_dir)
    print(f"[scycdb] {len(inputs)} panel genomes; {len(id2gene)} id->gene entries", file=sys.stderr)

    rows = []
    for g, src in inputs:
        hits = args.outdir / f"{g}.scyc.tsv"
        if not hits.exists():
            mode = "blastx" if src.suffix.lower() == ".fna" else "blastp"
            cmd = ["diamond", mode, "-k", "1", "-e", EVALUE, "-p", THREADS,
                   "-d", str(DB), "-q", str(src), "-o", str(hits), "--quiet"]
            subprocess.run(cmd, check=True)
        fam_count: dict[str, int] = defaultdict(int)
        seen_q: set[str] = set()
        with open(hits) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 2:
                    continue
                q, subj = c[0], c[1]
                if q in seen_q:        # -k 1 keeps best hit; guard anyway
                    continue
                seen_q.add(q)
                fam = id2gene.get(subj)
                if fam:
                    fam_count[fam] += 1
        for fam, n in sorted(fam_count.items()):
            rows.append({"sample": g, "family": fam, "count": n})
        print(f"[scycdb] {g}: {sum(fam_count.values())} mapped hits, {len(fam_count)} families",
              file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "family", "count"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[scycdb] wrote {args.out} ({len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
