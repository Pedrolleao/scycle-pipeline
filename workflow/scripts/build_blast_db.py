#!/usr/bin/env python3
"""
build_blast_db.py — fetch reviewed UniProt sequences for every target with
blast_fallback: true in targets.yaml, split the BLAST-gated targets into
their own database, and run makeblastdb (which also serves DIAMOND blastp
via the same FASTAs).

Two output FASTAs:
  resources/blast_db/unstable_refs.fasta    (Pfam-but-blast_fallback targets)
  resources/blast_db/blast_gated_refs.fasta (no-Pfam OR requires_blast_for_confirmation)

Header convention mirrors ESP_Search:  >{target_id}||{uniprot_accession}
This lets downstream code split on `||` to map a hit back to a target.

Usage:  python workflow/scripts/build_blast_db.py [--force]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from _common import (
    ROOT, load_targets, blast_fallback_targets, blast_gated_targets,
)

BLAST_DIR = ROOT / "resources" / "blast_db"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


def fetch_uniprot_fasta(accessions: list[str]) -> str:
    """Fetch the FASTA for a list of UniProt accessions in one request.

    Silently filters out UniParc cluster IDs (UPI*) — they're not fetchable
    via /uniprotkb/, and a single UPI in the batch causes HTTP 400 for the
    whole request. Curators sometimes paste them in from UniRef expansion.
    """
    if not accessions:
        return ""
    upi = [a for a in accessions if a.startswith("UPI")]
    accessions = [a for a in accessions if not a.startswith("UPI")]
    if upi:
        print(f"  ! filtered {len(upi)} UniParc ID(s) (UPI*) from UniProt batch: "
              f"{', '.join(upi[:5])}{'...' if len(upi) > 5 else ''}",
              file=sys.stderr)
    if not accessions:
        return ""
    query = " OR ".join(f"accession:{a}" for a in accessions)
    params = urllib.parse.urlencode({
        "query": query,
        "format": "fasta",
        "size": "500",
    })
    url = f"{UNIPROT_BASE}/search?{params}"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return r.read().decode()
    except Exception as e:
        print(f"  ! UniProt fetch failed: {e}", file=sys.stderr)
        return ""


def fetch_uniprot_by_gene(gene: str, organism: str | None = None,
                          limit: int = 5) -> str:
    """Pull up to `limit` reviewed entries by gene name (+ organism)."""
    terms = [f"gene:{gene}", "reviewed:true"]
    if organism:
        terms.append(f"organism_name:\"{organism}\"")
    params = urllib.parse.urlencode({
        "query": " AND ".join(terms),
        "format": "fasta",
        "size": str(limit),
    })
    url = f"{UNIPROT_BASE}/search?{params}"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return r.read().decode()
    except Exception as e:
        print(f"  ! UniProt gene fetch failed for {gene}: {e}", file=sys.stderr)
        return ""


def retag_headers(fasta: str, target_id: str) -> str:
    """Prepend `target_id||` to every UniProt header so downstream code can
    trace hits back to the target."""
    out = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            # UniProt headers look like ">sp|Q9X4K0|XX_YYY …"
            parts = line[1:].split("|", 2)
            acc = parts[1] if len(parts) >= 2 else parts[0]
            out.append(f">{target_id}||{acc} {line[1:]}")
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def build_one(out_fasta: Path, targets: list[dict]) -> None:
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    chunks = []
    for t in targets:
        accs = t.get("blast_refs_uniprot") or []
        if not accs:
            # AUDIT FIX (2026-05-27): do NOT fall back to a gene-NAME UniProt fetch. That pulled
            # cross-kingdom homonyms into the live DB (qmoA→Macaca SOD1/RNF146; sor→C. elegans
            # Sop-2 / maize Derlin; sorA→Penicillium polyketide synthase; phsA→Streptomyces
            # O-aminophenol oxidase). Targets with empty blast_refs_uniprot contribute NO seeds and
            # rely on KO/Pfam/custom-HMM + the S6 "BLAST-only is non-calling" rule. To enable BLAST
            # gating for such a target, CURATE function-verified accessions (the fetch_uniprot_by_gene
            # helper is retained but intentionally unused). NB: requires_blast_for_confirmation targets
            # MUST have curated seeds or they will always disqualify — currently only phsA/ttrA, now curated.
            print(f"  {t['id']}: no curated blast_refs_uniprot — skipping (no name-fetch).", flush=True)
            continue
        print(f"  fetching {t['id']}: {len(accs)} accession(s)", flush=True)
        fasta = fetch_uniprot_fasta(accs)
        if not fasta.strip():
            print(f"    ! no sequences returned for {t['id']}", file=sys.stderr)
            continue
        chunks.append(retag_headers(fasta, t["id"]))
        time.sleep(0.5)   # be polite to UniProt
    out_fasta.write_text("".join(chunks))
    if not out_fasta.stat().st_size:
        print(f"  ! {out_fasta} is empty — skipping makeblastdb", file=sys.stderr)
        return
    subprocess.run(["makeblastdb", "-in", str(out_fasta),
                    "-dbtype", "prot",
                    "-out", str(out_fasta.with_suffix("")),
                    "-title", out_fasta.stem],
                   check=True)
    # DIAMOND db (for the actual blast step in protein_mode.smk).
    subprocess.run(["diamond", "makedb", "--in", str(out_fasta),
                    "--db", str(out_fasta.with_suffix(".dmnd")),
                    "--quiet"],
                   check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    targets = load_targets()["targets"]
    unstable = [t for t in blast_fallback_targets(targets)
                if (t.get("pfam") or [])]
    gated = blast_gated_targets(targets)

    unstable_fa = BLAST_DIR / "unstable_refs.fasta"
    gated_fa = BLAST_DIR / "blast_gated_refs.fasta"

    if not args.force and unstable_fa.with_suffix(".phr").exists() \
            and gated_fa.with_suffix(".phr").exists():
        print("[build_blast_db] BLAST DBs already present; pass --force to rebuild")
        return

    print(f"[build_blast_db] unstable_refs.fasta: {len(unstable)} target(s)")
    build_one(unstable_fa, unstable)
    print(f"[build_blast_db] blast_gated_refs.fasta: {len(gated)} target(s)")
    build_one(gated_fa, gated)
    print("[build_blast_db] done.")


if __name__ == "__main__":
    main()
