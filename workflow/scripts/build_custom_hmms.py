#!/usr/bin/env python3
"""
build_custom_hmms.py — seed-list + automated-expansion HMM builder for
Tier 2 / Tier 3 targets in the e-waste pipeline.

Two sub-commands keep the workflow linear with a manual curation step in
between:

  expand --target {id}
    Read targets/{id}/manifest.yaml -> seed_accessions
    Fetch seed sequences from UniProt -> seeds.fasta
    Run phmmer against the cached SwissProt FASTA
    Write all hits to expanded.fasta (header retagged to {target_id}||{acc})
    Append each hit to manifest.yaml > expanded_candidates with keep: null
    --> Curator now opens the manifest and sets keep: true on rows to keep.

  build --target {id}
    Read manifest.yaml -> keep the rows with keep: true
    Materialize refs.fasta (subset of expanded.fasta)
    CD-HIT cluster at 0.7 identity (one rep per cluster)
    MAFFT --auto -> refs.aln
    hmmbuild --name {target_id} -> {target_id}.hmm
    Update manifest.yaml > n_seqs_final + curation_date

Reuses the UniProt fetch + header-retag pattern from build_blast_db.py.

Usage:
  python workflow/scripts/build_custom_hmms.py expand --target cyc2
  python workflow/scripts/build_custom_hmms.py build  --target cyc2
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import yaml

from _common import ROOT
from build_blast_db import fetch_uniprot_fasta, retag_headers

TARGETS_DIR = ROOT / "targets"
# UniRef90 (clustered at 90% identity) is the global search space for
# expand. SwissProt was too narrow for environmentally-relevant proteins
# like Cyc2 (TrEMBL-only family with no clean Pfam).
UNIREF90_GZ = ROOT / "resources" / ".cache" / "uniref90.fasta.gz"
UNIREF90_CACHE = ROOT / "resources" / ".cache" / "uniref90.fasta"
UNIREF90_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
)

MIN_SEQS = 8   # floor from WORK_PLAN section 4


def _target_dir(target_id: str) -> Path:
    d = TARGETS_DIR / target_id
    if not d.exists():
        sys.exit(f"error: {d} does not exist — create it and seed manifest.yaml first")
    return d


def _load_manifest(target_dir: Path) -> dict:
    mf = target_dir / "manifest.yaml"
    if not mf.exists():
        sys.exit(f"error: {mf} missing")
    with open(mf) as fh:
        return yaml.safe_load(fh) or {}


def _save_manifest(target_dir: Path, data: dict) -> None:
    mf = target_dir / "manifest.yaml"
    with open(mf, "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False)


# ── UniRef90 cache ──────────────────────────────────────────────────────────

def ensure_uniref90() -> Path:
    """Download UniRef90 once and cache it.

    Mirrors the Pfam-A pattern in build_hmm_db.py:59-72: keep the .gz on
    disk so re-extracting is free, and decompress to a sibling .fasta.
    UniRef90 is ~25 GB compressed, ~75 GB uncompressed — a one-off cost.
    """
    UNIREF90_CACHE.parent.mkdir(parents=True, exist_ok=True)
    if UNIREF90_CACHE.exists() and UNIREF90_CACHE.stat().st_size > 1_000_000_000:
        return UNIREF90_CACHE
    if not UNIREF90_GZ.exists():
        print(f"[build_custom_hmms] downloading UniRef90 → {UNIREF90_GZ} "
              "(one-off; ~47 GB compressed, ~150 GB uncompressed)…", flush=True)
        urllib.request.urlretrieve(UNIREF90_URL, UNIREF90_GZ)
    print(f"[build_custom_hmms] decompressing → {UNIREF90_CACHE}…", flush=True)
    with gzip.open(UNIREF90_GZ, "rb") as src, open(UNIREF90_CACHE, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    return UNIREF90_CACHE


# ── parsing helpers ─────────────────────────────────────────────────────────

def _bare_acc(tname: str) -> str:
    """Normalise a FASTA target name to a bare accession.

    Accepts SwissProt (`sp|Q9X4K0|XX_YYY`), TrEMBL (`tr|...|...`),
    UniRef (`UniRef90_A0A123`), retagged (`tid||acc`), or plain accessions.
    """
    if "||" in tname:
        return tname.split("||", 1)[1].split()[0]
    if "|" in tname:
        parts = tname.split("|")
        return parts[1] if len(parts) >= 2 else parts[0]
    if tname.startswith(("UniRef90_", "UniRef100_", "UniRef50_")):
        return tname.split("_", 1)[1]
    return tname.split()[0]


def _parse_fasta_str(text: str) -> dict[str, tuple[str, str]]:
    """Return {bare_acc: (header_line, sequence)} for a FASTA string."""
    out: dict[str, tuple[str, str]] = {}
    acc = None
    header = ""
    seq_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                out[acc] = (header, "".join(seq_lines))
            header = line.rstrip("\n")
            acc = _bare_acc(header[1:].split(None, 1)[0])
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if acc is not None:
        out[acc] = (header, "".join(seq_lines))
    return out


def _parse_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Return {bare_acc: (header_line, sequence)} for a FASTA file."""
    with open(path) as fh:
        return _parse_fasta_str(fh.read())


def _extract_records_by_acc(
    path: Path, wanted: set[str]
) -> dict[str, tuple[str, str]]:
    # Single-pass streaming filter — required for UniRef90 (~84 GB
    # uncompressed). _parse_fasta(uniref90) on a 30 GB-RAM box OOMs;
    # this keeps only records whose accession is in `wanted`.
    out: dict[str, tuple[str, str]] = {}
    remaining = set(wanted)
    keep = False
    acc = None
    header = ""
    seq_lines: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and acc is not None:
                    out[acc] = (header, "".join(seq_lines))
                    remaining.discard(acc)
                    if not remaining:
                        return out
                header = line.rstrip("\n")
                acc = _bare_acc(header[1:].split(None, 1)[0])
                keep = acc in remaining
                seq_lines = []
            elif keep:
                seq_lines.append(line.strip())
        if keep and acc is not None:
            out[acc] = (header, "".join(seq_lines))
    return out


def _organism_from_header(header: str) -> str:
    """Extract organism from either a UniProt (`OS=…`) or UniRef (`Tax=…`) header."""
    if "OS=" in header:
        rest = header.split("OS=", 1)[1]
        end = len(rest)
        for tag in (" OX=", " GN=", " PE=", " SV="):
            i = rest.find(tag)
            if i != -1 and i < end:
                end = i
        return rest[:end].strip()
    if "Tax=" in header:
        rest = header.split("Tax=", 1)[1]
        end = len(rest)
        for tag in (" TaxID=", " RepID=", " n="):
            i = rest.find(tag)
            if i != -1 and i < end:
                end = i
        return rest[:end].strip()
    return ""


def _parse_phmmer_tblout(path: Path) -> list[tuple[str, float]]:
    """Yield (bare_acc, full_bitscore) from a phmmer --tblout."""
    rows: list[tuple[str, float]] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                score = float(parts[5])
            except ValueError:
                continue
            rows.append((_bare_acc(parts[0]), score))
    return rows


# ── expand sub-command ──────────────────────────────────────────────────────

def cmd_expand(args: argparse.Namespace) -> None:
    target = args.target
    tdir = _target_dir(target)
    manifest = _load_manifest(tdir)
    seeds = manifest.get("seed_accessions") or []
    seeds = [s for s in seeds if s and not s.startswith("ACCESSION")]
    if not seeds:
        sys.exit(f"error: seed_accessions empty in {tdir/'manifest.yaml'}")
    print(f"[expand] {target}: fetching {len(seeds)} seed sequence(s) from UniProt")
    seeds_fa = fetch_uniprot_fasta(seeds)
    if not seeds_fa.strip():
        sys.exit(f"error: UniProt returned 0 sequences for seeds {seeds}")
    (tdir / "seeds.fasta").write_text(seeds_fa)

    sp = ensure_uniref90()

    tbl = tdir / ".phmmer.tblout"
    if args.reuse_tblout and tbl.exists():
        print(f"[expand] {target}: --reuse-tblout set and {tbl.name} exists; "
              "skipping phmmer")
    else:
        print(f"[expand] {target}: running phmmer vs UniRef90 (T={args.bitscore_threshold})")
        subprocess.run(
            ["phmmer",
             "--cpu", str(args.threads),
             "--noali",
             "-T", str(args.bitscore_threshold),
             "--tblout", str(tbl),
             str(tdir / "seeds.fasta"),
             str(sp)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    hits = _parse_phmmer_tblout(tbl)
    # Deduplicate: keep the best bitscore per accession
    best: dict[str, float] = {}
    for acc, score in hits:
        if acc not in best or score > best[acc]:
            best[acc] = score
    # Skip seeds — they are already trusted positives
    for s in seeds:
        best.pop(s, None)
    print(f"[expand] {target}: {len(best)} unique non-seed phmmer hit(s) above T={args.bitscore_threshold}")
    if not best:
        sys.exit("error: phmmer produced no hits above threshold — lower -T or revise seeds")

    # Pull full FASTAs for the hits (so the curator can grep organism names
    # straight out of expanded.fasta if they want to).
    hit_accs = sorted(best.keys(), key=lambda a: -best[a])
    if args.max_candidates:
        hit_accs = hit_accs[: args.max_candidates]
    expanded_fa = ""
    sp_records = _extract_records_by_acc(sp, set(hit_accs))
    for acc in hit_accs:
        rec = sp_records.get(acc)
        if rec is None:
            continue
        expanded_fa += f">{target}||{acc} {rec[0][1:]}\n"
        # wrap sequence at 60 cols for readability
        seq = rec[1]
        expanded_fa += "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
    # Also include seeds at the top of expanded.fasta so the curator sees them
    # as kept-by-default rows.
    seed_block = retag_headers(seeds_fa, target)
    (tdir / "expanded.fasta").write_text(seed_block + expanded_fa)

    # Rewrite manifest.expanded_candidates: seeds first (keep: true), then hits.
    cand_rows: list[dict] = []
    seed_records = _parse_fasta(tdir / "seeds.fasta")
    for s in seeds:
        rec = seed_records.get(s)
        if rec is None:
            continue
        cand_rows.append({
            "acc": s,
            "organism": _organism_from_header(rec[0]),
            "length": len(rec[1]),
            "phmmer_score": None,
            "keep": True,
            "note": "seed",
        })
    for acc in hit_accs:
        rec = sp_records.get(acc)
        if rec is None:
            continue
        cand_rows.append({
            "acc": acc,
            "organism": _organism_from_header(rec[0]),
            "length": len(rec[1]),
            "phmmer_score": round(best[acc], 1),
            "keep": None,
            "note": "",
        })
    manifest["expanded_candidates"] = cand_rows
    _save_manifest(tdir, manifest)
    tbl.unlink(missing_ok=True)
    print(f"[expand] {target}: wrote expanded.fasta and {len(cand_rows)} "
          f"candidates to manifest.yaml — review and set keep: true on the ones to keep")


# ── build sub-command ───────────────────────────────────────────────────────

def cmd_build(args: argparse.Namespace) -> None:
    target = args.target
    tdir = _target_dir(target)
    manifest = _load_manifest(tdir)
    cands = manifest.get("expanded_candidates") or []
    kept = [c for c in cands if c.get("keep") is True]
    # MIN_SEQS default 8, override via manifest.min_seqs_override when the
    # family is biologically narrow (rare positives + extensive negatives is
    # still a working HMM). The override must be documented in manifest.notes.
    min_floor = int(manifest.get("min_seqs_override") or MIN_SEQS)
    if len(kept) < min_floor:
        sys.exit(f"error: only {len(kept)} sequence(s) marked keep: true "
                 f"(floor is {min_floor}). Either curate more rows or move "
                 f"this target to Tier 3 (BLAST-only).")
    print(f"[build] {target}: {len(kept)} curated sequence(s) → CD-HIT → MAFFT → hmmbuild")

    # Materialize refs.fasta from expanded.fasta, subset to kept accs.
    # Curators occasionally add accessions they found via fresh UniProt query
    # rather than from the expand step's output (happens when the seed was
    # wrong and the expansion was poisoned — e.g. hcnA's Q9I3F9 being an
    # alpha/beta hydrolase, not hcnA). Auto-fetch missing accessions from
    # UniProt so the manifest is the source of truth, not expanded.fasta.
    exp_records = _parse_fasta(tdir / "expanded.fasta")
    missing_accs = [c["acc"] for c in kept if c["acc"] not in exp_records]
    if missing_accs:
        print(f"[build] {target}: {len(missing_accs)} kept accession(s) not in "
              f"expanded.fasta — fetching from UniProt", file=sys.stderr)
        fresh = fetch_uniprot_fasta(missing_accs)
        if fresh.strip():
            fresh_records = _parse_fasta_str(fresh)
            exp_records.update(fresh_records)
    refs_chunks = []
    for c in kept:
        rec = exp_records.get(c["acc"])
        if rec is None:
            print(f"  ! kept accession {c['acc']} not in expanded.fasta and "
                  f"not fetchable from UniProt — skipping", file=sys.stderr)
            continue
        refs_chunks.append(f">{target}||{c['acc']}\n{rec[1]}\n")
    refs_fa = tdir / "refs.fasta"
    refs_fa.write_text("".join(refs_chunks))

    # CD-HIT redundancy clustering. Default 0.7 per WORK_PLAN §4, but highly
    # conserved narrow families (e.g. rusticyanin) collapse below the MIN_SEQS
    # floor at 0.7 — raise per-target via manifest.cdhit_threshold.
    cdhit_c = float(manifest.get("cdhit_threshold") or 0.7)
    cdhit_n = 5 if cdhit_c >= 0.7 else (4 if cdhit_c >= 0.6 else 3)
    cdhit_out = tdir / "refs.cdhit.fasta"
    subprocess.run(
        ["cd-hit", "-i", str(refs_fa), "-o", str(cdhit_out),
         "-c", f"{cdhit_c}", "-n", str(cdhit_n),
         "-T", str(args.threads), "-M", "0"],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    # Replace refs.fasta with the clustered version.
    refs_fa.write_text(cdhit_out.read_text())
    cdhit_out.unlink(missing_ok=True)
    (tdir / "refs.fasta.clstr").unlink(missing_ok=True)
    Path(str(cdhit_out) + ".clstr").unlink(missing_ok=True)
    n_after = sum(1 for line in refs_fa.read_text().splitlines() if line.startswith(">"))
    if n_after < min_floor:
        sys.exit(f"error: only {n_after} sequence(s) after CD-HIT clustering "
                 f"(floor is {min_floor}). Curate more diverse sequences "
                 f"or raise cdhit_threshold toward 0.95.")

    # MAFFT alignment. --anysymbol handles selenocysteine (U) and pyrrolysine
    # (O), which appear in legitimate sequences like fdhF (E. coli formate
    # dehydrogenase F has Sec at residue 140).
    aln = tdir / "refs.aln"
    print(f"[build] {target}: aligning {n_after} sequence(s) with MAFFT")
    with open(aln, "w") as out:
        subprocess.run(
            ["mafft", "--anysymbol", "--auto", "--thread", str(args.threads), str(refs_fa)],
            check=True,
            stdout=out,
            stderr=subprocess.DEVNULL,
        )

    # hmmbuild.
    hmm = tdir / f"{target}.hmm"
    subprocess.run(
        ["hmmbuild", "--cpu", str(args.threads), "-n", target,
         str(hmm), str(aln)],
        check=True,
        stdout=subprocess.DEVNULL,
    )

    manifest["n_seqs_final"] = n_after
    manifest["curation_date"] = _dt.date.today().isoformat()
    _save_manifest(tdir, manifest)
    print(f"[build] {target}: wrote {hmm} ({n_after} sequence(s))")


# ── main ────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_exp = sub.add_parser("expand", help="phmmer-expand seeds against SwissProt")
    p_exp.add_argument("--target", required=True)
    p_exp.add_argument("--threads", type=int, default=4)
    p_exp.add_argument("--bitscore-threshold", type=float, default=80.0,
                       help="phmmer -T threshold (default 80; lower → broader)")
    p_exp.add_argument("--max-candidates", type=int, default=200,
                       help="cap the number of phmmer candidates written back "
                            "(default 200; null/0 = unlimited)")
    p_exp.add_argument("--reuse-tblout", action="store_true",
                       help="if targets/{id}/.phmmer.tblout already exists, "
                            "skip the phmmer run and parse the cached output. "
                            "Use for recovery after a post-phmmer crash.")
    p_exp.set_defaults(func=cmd_expand)

    p_build = sub.add_parser("build", help="CD-HIT + MAFFT + hmmbuild from kept rows")
    p_build.add_argument("--target", required=True)
    p_build.add_argument("--threads", type=int, default=4)
    p_build.set_defaults(func=cmd_build)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
