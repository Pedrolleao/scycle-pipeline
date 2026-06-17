#!/usr/bin/env python3
"""
calibrate_tc.py — set the trusted-cutoff (TC) bitscore for a custom HMM.

Inputs:
  targets/{target_id}/refs.fasta      (positives — same set used to build the HMM)
  targets/{target_id}/{target_id}.hmm (the built HMM)
  targets/{target_id}/manifest.yaml > negatives (UniProt accessions)

Outputs:
  targets/{target_id}/negatives.fasta
  targets/{target_id}/tc_calibration.tsv  (label, acc, bitscore, evalue)
  manifest.yaml > tc_bitscore + tc_rationale

Algorithm:
  1. hmmsearch the HMM against refs.fasta → positive bitscore distribution
  2. hmmsearch the HMM against negatives.fasta → negative bitscore distribution
  3. Pick TC at the *largest gap* between the lowest positive and highest
     negative score: TC = (min_positive + max_negative) / 2 if min_pos > max_neg,
     otherwise warn and fall back to min(positives) * 0.95.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

from _common import ROOT
from build_blast_db import fetch_uniprot_fasta, retag_headers

TARGETS_DIR = ROOT / "targets"


def _load_manifest(tdir: Path) -> dict:
    with open(tdir / "manifest.yaml") as fh:
        return yaml.safe_load(fh) or {}


def _save_manifest(tdir: Path, data: dict) -> None:
    with open(tdir / "manifest.yaml", "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False)


def _hmmsearch_bitscores(hmm: Path, fasta: Path) -> dict[str, float]:
    """Return {target_name: best_full_bitscore} from hmmsearch --tblout."""
    if not fasta.exists() or fasta.stat().st_size == 0:
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".tblout", delete=False) as tf:
        tbl = Path(tf.name)
    subprocess.run(
        ["hmmsearch", "--noali", "--tblout", str(tbl), str(hmm), str(fasta)],
        check=True, stdout=subprocess.DEVNULL,
    )
    out: dict[str, float] = {}
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            tname = parts[0]
            try:
                score = float(parts[5])
            except ValueError:
                continue
            if tname not in out or score > out[tname]:
                out[tname] = score
    tbl.unlink(missing_ok=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--margin-fraction", type=float, default=0.05,
                    help="Fallback margin when positives and negatives overlap "
                         "(default 0.05 → TC = min_positive * 0.95)")
    args = ap.parse_args()

    tdir = TARGETS_DIR / args.target
    if not tdir.exists():
        sys.exit(f"error: {tdir} does not exist")
    manifest = _load_manifest(tdir)

    hmm = tdir / f"{args.target}.hmm"
    refs = tdir / "refs.fasta"
    if not hmm.exists() or not refs.exists():
        sys.exit(f"error: missing {hmm} or {refs} — run build_custom_hmms.py build first")

    # Fetch negatives.
    negs = manifest.get("negatives") or []
    negs = [n for n in negs if n and not n.startswith("ACCESSION")]
    if not negs:
        sys.exit(f"error: negatives empty in manifest.yaml — add ~15 close-but-not-{args.target} accessions")
    print(f"[calibrate_tc] {args.target}: fetching {len(negs)} negative(s) from UniProt")
    neg_fa = fetch_uniprot_fasta(negs)
    if not neg_fa.strip():
        sys.exit(f"error: UniProt returned 0 sequences for negatives {negs}")
    (tdir / "negatives.fasta").write_text(
        retag_headers(neg_fa, f"NEG_{args.target}")
    )

    pos_scores = _hmmsearch_bitscores(hmm, refs)
    neg_scores = _hmmsearch_bitscores(hmm, tdir / "negatives.fasta")
    if not pos_scores:
        sys.exit("error: HMM scored 0 of its own positives — refs.fasta corrupt?")

    print(f"[calibrate_tc] {args.target}: {len(pos_scores)} positives "
          f"(min={min(pos_scores.values()):.1f}, max={max(pos_scores.values()):.1f})")
    if neg_scores:
        print(f"[calibrate_tc] {args.target}: {len(neg_scores)} negatives "
              f"(min={min(neg_scores.values()):.1f}, max={max(neg_scores.values()):.1f})")
    else:
        print(f"[calibrate_tc] {args.target}: no negatives scored above default e-value "
              "(this is the desired outcome — HMM is highly specific)")

    min_pos = min(pos_scores.values())
    max_neg = max(neg_scores.values()) if neg_scores else 0.0
    if min_pos > max_neg:
        tc = (min_pos + max_neg) / 2.0
        rationale = (f"clean gap: min_positive={min_pos:.1f} > max_negative={max_neg:.1f}; "
                     f"TC = midpoint = {tc:.1f}")
    else:
        tc = round(min_pos * (1.0 - args.margin_fraction), 1)
        rationale = (f"OVERLAP: min_positive={min_pos:.1f} <= max_negative={max_neg:.1f}; "
                     f"TC = min_positive * (1 - {args.margin_fraction}) = {tc:.1f}. "
                     f"WARNING: HMM is not fully discriminating — widen negatives or curate more positives.")
        print(f"[calibrate_tc] WARNING: positives and negatives overlap", file=sys.stderr)

    # Write tc_calibration.tsv.
    with open(tdir / "tc_calibration.tsv", "w") as fh:
        fh.write("label\tname\tbitscore\n")
        for name, score in sorted(pos_scores.items(), key=lambda x: -x[1]):
            fh.write(f"positive\t{name}\t{score:.1f}\n")
        for name, score in sorted(neg_scores.items(), key=lambda x: -x[1]):
            fh.write(f"negative\t{name}\t{score:.1f}\n")

    manifest["tc_bitscore"] = round(float(tc), 1)
    manifest["tc_rationale"] = rationale
    _save_manifest(tdir, manifest)
    print(f"[calibrate_tc] {args.target}: TC = {tc:.1f} — written to manifest.yaml")
    print(f"  rationale: {rationale}")


if __name__ == "__main__":
    main()
