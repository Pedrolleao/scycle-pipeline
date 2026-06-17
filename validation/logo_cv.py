#!/usr/bin/env python3
"""
logo_cv.py — leave-one-clade-out cross-validation for a custom HMM (SP3).

Question: does a lifted custom HMM detect its target by capturing the conserved fold
(generalizes to clades NOT in its training refs), or is detection circular (it only
works because the test genome's homolog was a training ref)?

Method: read targets/{id}/manifest.yaml kept refs + their organisms, drop every ref whose
organism matches a held-out clade keyword, rebuild the HMM (MAFFT --auto -> hmmbuild) from
the remaining refs, then hmmsearch the rebuilt HMM against each test proteome and report the
best full-sequence bitscore vs the re-calibrated TC. A held-out positive that still scores
well above TC (and above the negatives) is genuine generalization.

Note: the SP1 panel already provides an in-vivo leave-CLASS-out test — the independent SOB
positives (Sulfurimonas/epsilon, Thiobacillus/beta) belong to proteobacterial classes ABSENT
from the soxB/soxD training refs (alpha + gamma + Chlorobi + Thermus), yet are detected. This
script makes that rigorous and extends it to the represented (alpha) clade.

Usage:
  python validation/logo_cv.py --target soxB \
      --holdout Paracoccus Ancylobacter Gemmobacter Cereibacter Rhodovulum Ruegeria Paracoccaceae \
      --test Pdenitrificans_PD1222 Dshibae_DFL12 Rpalustris_CGA009 \
             Tdenitrificans_ATCC25259 Sdenitrificans_DSM1251 Bsubtilis_168 Smeliloti_1021
  (proteomes resolved from --panel, default ../sp1_panel; needs mafft + hmmbuild + hmmsearch)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def kept_refs(target: str) -> dict[str, str]:
    mf = yaml.safe_load(open(ROOT / "targets" / target / "manifest.yaml"))
    return {c["acc"]: c.get("organism", "")
            for c in mf.get("expanded_candidates", []) if c.get("keep") is True}


def tc_of(target: str) -> float:
    for line in (ROOT / "resources" / "hmm" / "tc_cutoffs.tsv").read_text().splitlines():
        p = line.split("\t")
        if p and p[0] == target and len(p) >= 3 and p[2]:
            return float(p[2])
    return float("nan")


def best_score(hmm: Path, faa: Path, work: Path) -> float | None:
    tbl = work / "out.tbl"
    subprocess.run(["hmmsearch", "--noali", "--tblout", str(tbl), str(hmm), str(faa)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    if not tbl.exists():
        return None
    for line in tbl.read_text().splitlines():
        if not line.startswith("#"):
            return float(line.split()[5])   # full-seq score column
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--holdout", nargs="+", required=True,
                    help="organism-name keywords; refs whose organism matches any are dropped")
    ap.add_argument("--test", nargs="+", required=True, help="panel proteome stems to score")
    ap.add_argument("--panel", type=Path, default=ROOT.parent / "sp1_panel")
    args = ap.parse_args()

    refs = kept_refs(args.target)
    held = {a for a, o in refs.items() if any(k in o for k in args.holdout)}
    print(f"[logo] {args.target}: dropping {len(held)}/{len(refs)} refs matching {args.holdout}",
          file=sys.stderr)

    # subset refs.fasta to the non-held-out refs
    seqs, acc = {}, None
    for line in (ROOT / "targets" / args.target / "refs.fasta").read_text().splitlines(keepends=True):
        if line.startswith(">"):
            acc = line.split("||")[-1].strip()
            seqs[acc] = [line]
        elif acc:
            seqs[acc].append(line)
    kept = [a for a in seqs if a not in held]

    tc = tc_of(args.target)
    with tempfile.TemporaryDirectory() as wd:
        work = Path(wd)
        (work / "na.fasta").write_text("".join(l for a in kept for l in seqs[a]))
        subprocess.run(f"mafft --auto {work}/na.fasta > {work}/na.aln",
                       shell=True, stderr=subprocess.DEVNULL, check=True)
        subprocess.run(["hmmbuild", "--amino", "-n", f"{args.target}_logo",
                        str(work / "logo.hmm"), str(work / "na.aln")],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print(f"\nleave-out {args.target} HMM ({len(kept)} refs) vs panel — TC={tc}")
        print(f"  {'proteome':<28}{'bitscore':>10}  call")
        for stem in args.test:
            faa = args.panel / f"{stem}.faa"
            s = best_score(work / "logo.hmm", faa, work) if faa.exists() else None
            call = "—" if s is None else (">=TC ✓" if s >= tc else "below TC")
            print(f"  {stem:<28}{(f'{s:.1f}' if s else 'none'):>10}  {call}")


if __name__ == "__main__":
    main()
