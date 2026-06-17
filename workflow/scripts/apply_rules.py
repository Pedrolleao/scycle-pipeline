#!/usr/bin/env python3
"""
apply_rules.py — combine hmmscan + DIAMOND-blastp evidence, emit one row per
target into calls/scycle_calls.tsv.

Adapted from /home/dmin/Research/Asgard_Vault/ESP_Search/pipeline/apply_rules.py
(evaluate_rule, lines 67–127). Key changes:
  - Input is hmmscan domtblout (replaces InterProScan TSV)
  - "Narrow cluster hit" is replaced by "BLAST hit against a curated UniProt
    reference tagged for this target" (header convention: >{target_id}||{acc})
  - Status set is unchanged: confirmed | domain-only | narrow-no-IPR | disqualified
  - CadA gets a special disqualification: PF00330 hit without a CadA BLAST hit
    is treated as aconitase (disqualified)

Two modes (--mode protein | --mode read):
  - protein: evidence from hmmscan + DIAMOND-blastp
  - read   : evidence from a presence TSV produced by parse_coverage.py
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

KO_RE = re.compile(r"^K\d{4,}$")


# ───────────────────────────── parsers ───────────────────────────────────────

def parse_hmmscan_domtbl(path: Path, tc_cutoffs: dict[str, float],
                         qcov_cutoffs: dict[str, float] | None = None
                        ) -> dict[str, dict]:
    """Return {protein_id: {pfam: {pfam_id: best_evalue},
                             custom: {target_id: best_evalue}}}.

    Pfam profiles are identified by target_acc starting with "PF" (e.g.
    PF13435.10). Custom HMMs (built via build_custom_hmms.py with
    `hmmbuild --name target_id`) carry the target_id in target_name and
    typically have no ACC field (`-`).

    `qcov_cutoffs` is a per-target {target_id: min_query_coverage} dict.
    When set, domains with env-coverage of the query protein below the
    threshold are dropped — prevents small-domain HMMs (e.g. 104-aa hcnA)
    from false-positive-hitting the matching domain of a large multidomain
    protein (~950 aa Fe-S-containing oxidoreductase).
    """
    qcov_cutoffs = qcov_cutoffs or {}
    out: dict[str, dict] = defaultdict(lambda: {"pfam": {}, "custom": {}, "ko": {}})
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            target_name = parts[0]
            target_acc = parts[1]                 # "PF#####.NN" or "-"
            qlen = int(parts[5])
            query_name = parts[3]
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            env_from = int(parts[19])
            env_to = int(parts[20])
            if target_acc.startswith("PF"):
                key = target_acc.split(".")[0]
                bucket = "pfam"
            elif KO_RE.match(target_name) or KO_RE.match(target_acc):
                # KOfam profile — name (or acc) is the KO number, e.g. K02588.
                key = target_name if KO_RE.match(target_name) else target_acc
                bucket = "ko"
            else:
                key = target_name
                bucket = "custom"
            tc = tc_cutoffs.get(key)
            # Apply TC cutoff if available, else trust the e-value filter that
            # hmmscan already applied via -E.
            if tc is not None and full_score < tc:
                continue
            # Apply per-target query-coverage filter (defends against
            # small-domain HMMs over-calling on multidomain proteins).
            min_qcov = qcov_cutoffs.get(key)
            if min_qcov is not None and qlen > 0:
                qcov = (env_to - env_from + 1) / qlen
                if qcov < min_qcov:
                    continue
            entry = out[query_name][bucket]
            prev = entry.get(key)
            if prev is None or full_evalue < prev:
                entry[key] = full_evalue
    return out


def parse_blast_tsv(path: Path) -> dict[str, dict[str, list[tuple[str, float, float]]]]:
    """Return {protein_id: {target_id: [(uniprot_acc, pident, evalue), …]}}.
    Subject header convention: target_id||uniprot_acc."""
    out: dict[str, dict] = defaultdict(lambda: defaultdict(list))
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qid = parts[0]
            sid = parts[1]
            if "||" not in sid:
                continue
            target_id, acc = sid.split("||", 1)
            pident = float(parts[2])
            evalue = float(parts[10])
            out[qid][target_id].append((acc.split(" ", 1)[0], pident, evalue))
    return out


def load_tc_cutoffs(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2]:
                try:
                    out[parts[0]] = float(parts[2])
                except ValueError:
                    pass
    return out


# ───────────────────────────── per-target evaluation ─────────────────────────

def evaluate_target(target: dict,
                    pfam_hits: dict[str, dict[str, float]],
                    custom_hits: dict[str, dict[str, float]],
                    ko_hits: dict[str, dict[str, float]],
                    blast_hits: dict[str, dict],
                    pident_default: float, qcov_default: float,
                    ) -> tuple[str, dict] | None:
    """Return (status, evidence_dict) for the BEST protein matching this target,
    or None if nothing matches.

    Status set:
      confirmed     — sequence signature met (Pfam OR custom HMM) AND BLAST hit
                       (or signature-only when no fallback is configured)
      domain-only   — signature met, no BLAST hit (only meaningful when fallback set)
      narrow-no-IPR — BLAST hit, signature not met
      disqualified  — requires_blast_for_confirmation set, signature met, BLAST missing
    """
    tid = target["id"]
    is_custom = bool(target.get("custom_hmm"))
    req_pfams = set(target.get("pfam") or [])
    target_kos = set(target.get("ko") or [])
    logic = (target.get("pfam_logic") or "any").lower()   # any | all
    fallback = bool(target.get("blast_fallback"))
    requires_blast_for_confirmation = bool(
        target.get("requires_blast_for_confirmation"))
    pident_min = float(target.get("blast_identity_min")
                       or pident_default)

    # Collect per-protein evidence summaries.
    candidates = []
    proteins = (set(pfam_hits.keys()) | set(custom_hits.keys())
                | set(ko_hits.keys()) | set(blast_hits.keys()))
    for prot in proteins:
        p_hits = pfam_hits.get(prot, {})
        c_hits = custom_hits.get(prot, {})
        k_hits = ko_hits.get(prot, {})
        b_hits_all = blast_hits.get(prot, {}).get(tid, [])
        b_hits = [(a, pi, ev) for (a, pi, ev) in b_hits_all
                  if pi >= pident_min]
        # Sequence-signature check. Precedence: custom HMM > KOfam KO > Pfam.
        # The custom-HMM tier is preferred when a clade model has been built
        # (hardening phase); until then KO is the primary signature, so a
        # `custom_hmm: true` target with no built model still resolves via KO.
        sig_source = None
        sig_evalue = None
        present_pfams: set[str] = set()
        present_kos: set[str] = set()
        if is_custom and tid in c_hits:
            sig_source = "custom-hmm"
            sig_evalue = c_hits[tid]
        if sig_source is None and target_kos:
            present_kos = target_kos & set(k_hits.keys())
            if present_kos:
                sig_source = "ko"
                sig_evalue = min(k_hits[k] for k in present_kos)
        if sig_source is None and req_pfams:
            present_pfams = req_pfams & set(p_hits.keys())
            ok = ((logic == "all" and present_pfams == req_pfams)
                  or (logic == "any" and bool(present_pfams)))
            if ok:
                sig_source = "pfam"
                sig_evalue = min(p_hits[p] for p in present_pfams)
        sig_ok = sig_source is not None
        blast_ok = bool(b_hits)
        if not (sig_ok or blast_ok):
            continue
        # Status.
        if requires_blast_for_confirmation and sig_ok and not blast_ok:
            status = "disqualified"
        elif sig_ok and (fallback and blast_ok or not fallback):
            status = "confirmed"
        elif sig_ok and fallback and not blast_ok:
            status = "domain-only"
        elif not sig_ok and blast_ok:
            # BLAST-only hit with NO sequence signature. For any target that has a
            # primary signature tier (custom HMM, KOfam KO, or Pfam), a BLAST-only hit
            # means that signature was rejected → cross-reactivity to a generic
            # homologue (often a UniProt entry auto-fetched by a shared gene name),
            # not a real call (S6, 2026-05-25 hardening). Only true Tier-3 targets
            # (no custom HMM, no KO, no Pfam — e.g. otr) legitimately call on BLAST
            # alone, where narrow-no-IPR remains a meaningful "BLAST says yes" status.
            if is_custom or target_kos or req_pfams:
                continue
            status = "narrow-no-IPR"
        else:
            continue
        evidence_source = sig_source if sig_ok else "blast"
        best_b = min(b_hits, key=lambda x: x[2]) if b_hits else None
        candidates.append({
            "protein_id": prot,
            "status": status,
            "evidence_source": evidence_source,
            "pfam_hits": (sorted(present_kos) if sig_source == "ko"
                          else sorted(present_pfams) if sig_source == "pfam"
                          else [tid] if sig_source == "custom-hmm"
                          else []),
            "pfam_evalue": sig_evalue,
            "blast_acc": best_b[0] if best_b else "",
            "blast_pident": best_b[1] if best_b else "",
            "blast_evalue": best_b[2] if best_b else "",
        })
    if not candidates:
        return None
    # Prefer confirmed > domain-only > narrow-no-IPR > disqualified.
    rank = {"confirmed": 0, "domain-only": 1,
            "narrow-no-IPR": 2, "disqualified": 3}
    best = sorted(candidates,
                  key=lambda c: (rank[c["status"]],
                                 c["pfam_evalue"] or c["blast_evalue"] or 1e9))[0]
    return best["status"], best


# ───────────────────────────── read-mode shortcut ────────────────────────────

def evaluate_target_from_reads(target: dict,
                               read_presence: dict[str, dict]) -> tuple[str, dict] | None:
    """For Mode B, presence is per-target boolean from parse_coverage.py."""
    tid = target["id"]
    pr = read_presence.get(tid)
    if pr is None:
        return None
    if not pr.get("present"):
        return None
    status = "confirmed"   # if breadth + depth + identity met, treat as confirmed
    return status, {
        "protein_id": "",
        "status": status,
        "pfam_hits": [],
        "pfam_evalue": "",
        "blast_acc": pr.get("source_acc", ""),
        "blast_pident": pr.get("identity", ""),
        "blast_evalue": "",
        "breadth_pct": pr.get("breadth_pct"),
        "depth_median": pr.get("depth_median"),
    }


def load_read_presence(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not path or not Path(path).exists():
        return out
    import csv
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tid = row["target_id"]
            out[tid] = {
                "present": row["present"] in ("True", "true", "1"),
                "source_acc": row.get("source_acc", ""),
                "identity": float(row["identity_pct"] or 0.0),
                "breadth_pct": float(row["breadth_pct"] or 0.0),
                "depth_median": float(row["depth_median"] or 0.0),
            }
    return out


# ───────────────────────────── dsr direction ─────────────────────────────────

def parse_prodigal_coords(faa_path: Path | None) -> dict[str, tuple[str, int, int]]:
    """Parse gene coordinates from a Prodigal .faa. Prodigal headers are
    `>{seqid}_{n} # {start} # {end} # {strand} # ID=...`; the downstream protein
    id is `{seqid}_{n}` and the contig is the seqid (id minus the trailing _{n}).
    Returns {protein_id: (contig, start, end)}, or {} for a pre-called proteome
    (no `#`-delimited coordinate header) — synteny resolution is then skipped."""
    out: dict[str, tuple[str, int, int]] = {}
    if not faa_path or not Path(faa_path).exists():
        return out
    with open(faa_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = [p.strip() for p in line[1:].split("#")]
            if len(parts) < 4:
                continue                      # not a Prodigal coordinate header
            protid = parts[0].split()[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            out[protid] = (protid.rsplit("_", 1)[0], start, end)
    return out


def _present(rows_by_id: dict[str, dict], tid: str) -> bool:
    r = rows_by_id.get(tid)
    return r is not None and r["status"] in (
        "confirmed", "domain-only", "narrow-no-IPR")


def resolve_dsr_direction(rows_by_id: dict[str, dict],
                          coords: dict[str, tuple[str, int, int]] | None = None,
                          window: int = 20000) -> str | None:
    """Call the metabolic DIRECTION of a dsrAB-carrying genome.

    dsrA/dsrB are the SAME genes in dissimilatory sulfate reducers (REDUCTIVE Dsr,
    SO3²⁻→H2S) and in sulfur oxidizers (REVERSE/oxidative rDSR, H2S/S⁰→SO3²⁻).
    Sequence alone cannot separate the two directions — the genomic companions do:
      reductive  — dsrD present (and usually qmoABC / dsrMK); Sox machinery absent.
                   dsrD sits in the canonical dsrABD operon, so when gene coordinates
                   are available a dsrD ORF syntenic with dsrA/dsrB is a strong call.
      oxidative  — dsrD absent, Sox core / sqr present (reverse-Dsr sulfur oxidizer).
    Returns 'reductive' | 'oxidative' | 'ambiguous', or None if dsrAB is absent.
    Presence-based (works on pre-called proteomes); refined by dsrD↔dsrAB synteny
    when `coords` is supplied (nucleotide/MAG input)."""
    if not (_present(rows_by_id, "dsrA") or _present(rows_by_id, "dsrB")):
        return None
    has_dsrD = _present(rows_by_id, "dsrD")
    has_qmo = any(_present(rows_by_id, t) for t in ("qmoA", "qmoB", "qmoC"))
    has_sox = any(_present(rows_by_id, t)
                  for t in ("soxB", "soxA", "soxY", "soxX", "sqr", "soxC"))

    # Synteny refinement: a dsrD ORF adjacent to dsrA/dsrB (canonical dsrABD operon)
    # is a strong reductive signal. Requires gene coordinates (nucleotide/MAG input).
    syntenic_dsrD = False
    if coords:
        dsrD_p = (rows_by_id.get("dsrD") or {}).get("protein_id", "")
        if dsrD_p in coords:
            cd, sd, ed = coords[dsrD_p]
            lo1, hi1 = min(sd, ed), max(sd, ed)
            for t in ("dsrA", "dsrB"):
                p = (rows_by_id.get(t) or {}).get("protein_id", "")
                if p not in coords:
                    continue
                cc, sc, ec = coords[p]
                if cc == cd and max(lo1, min(sc, ec)) - min(hi1, max(sc, ec)) <= window:
                    syntenic_dsrD = True
                    break

    # Priority (fixed 2026-05-25 SP2): dsrD is the REDUCTIVE-specific marker; Sox/sqr
    # machinery (with dsrD absent) marks an oxidative reverse-Dsr sulfur oxidizer; qmo is
    # only a weak reductive fallback (it occurs in BOTH directions, so it must NOT override
    # the Sox signal — that bug mislabeled Thiobacillus, a sulfur oxidizer, as reductive).
    if syntenic_dsrD or has_dsrD:
        return "reductive"
    if has_sox:
        return "oxidative"
    if has_qmo:
        return "reductive"
    return "ambiguous"


# ───────────────────────────── main ──────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", choices=["protein", "read"], required=True)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--hmm", type=Path)
    ap.add_argument("--blast-unstable", type=Path)
    ap.add_argument("--blast-gated", type=Path)
    ap.add_argument("--read-presence", type=Path,
                    help="TSV from parse_coverage.py (Mode B)")
    ap.add_argument("--tc-cutoffs", type=Path,
                    default=Path("resources/hmm/tc_cutoffs.tsv"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pident-default", type=float, default=30.0)
    ap.add_argument("--qcov-default", type=float, default=50.0)
    ap.add_argument("--gene-coords", type=Path,
                    help="proteome .faa; if it carries Prodigal coordinate headers "
                         "(nucleotide/MAG input), refine the dsrAB reductive-vs-"
                         "oxidative direction call by dsrD↔dsrAB operon synteny. "
                         "Falls back to presence-based for pre-called proteomes.")
    args = ap.parse_args()

    with open(args.targets) as fh:
        targets = yaml.safe_load(fh)["targets"]

    if args.mode == "protein":
        tc = load_tc_cutoffs(args.tc_cutoffs)
        # Per-target hmm_min_qcov override from targets.yaml — defends small
        # domain HMMs (e.g. 104-aa hcnA) from hitting domains of large
        # multidomain proteins.
        qcov_cutoffs = {t["id"]: float(t["hmm_min_qcov"])
                        for t in targets if t.get("hmm_min_qcov") is not None}
        hmm_parsed = parse_hmmscan_domtbl(args.hmm, tc, qcov_cutoffs) if args.hmm else {}
        blast_hits = {}
        for p in (args.blast_unstable, args.blast_gated):
            if p and p.exists():
                for prot, td in parse_blast_tsv(p).items():
                    blast_hits.setdefault(prot, {}).update(td)
        # Split per signature tier: pfam_hits[prot]={pfam_id: evalue};
        # custom_hits[prot]={tid: evalue}; ko_hits[prot]={KO: evalue}.
        pfam_hits_flat = {p: d["pfam"] for p, d in hmm_parsed.items()}
        custom_hits_flat = {p: d["custom"] for p, d in hmm_parsed.items()}
        ko_hits_flat = {p: d["ko"] for p, d in hmm_parsed.items()}

        rows = []
        for t in targets:
            res = evaluate_target(t, pfam_hits_flat, custom_hits_flat,
                                  ko_hits_flat, blast_hits,
                                  args.pident_default, args.qcov_default)
            if res is None:
                rows.append({
                    "target_id": t["id"], "name": t["name"],
                    "category": t["category"], "complex": t.get("complex", ""),
                    "status": "absent", "evidence_source": "",
                    "protein_id": "", "pfam_hits": "",
                    "best_pfam_evalue": "", "blast_acc": "",
                    "blast_pident": "", "blast_evalue": "",
                })
                continue
            status, ev = res
            rows.append({
                "target_id": t["id"], "name": t["name"],
                "category": t["category"], "complex": t.get("complex", ""),
                "status": status,
                "evidence_source": ev.get("evidence_source", ""),
                "protein_id": ev["protein_id"],
                "pfam_hits": ",".join(ev["pfam_hits"]),
                "best_pfam_evalue": ev["pfam_evalue"],
                "blast_acc": ev["blast_acc"],
                "blast_pident": ev["blast_pident"],
                "blast_evalue": ev["blast_evalue"],
            })

        # dsrAB direction call: reductive (dissimilatory sulfate reduction) vs
        # reverse/oxidative rDSR (sulfur oxidizer). dsrD/qmo presence — and, when
        # gene coordinates are present, dsrD↔dsrAB synteny — ⇒ reductive; dsrD absent
        # with Sox/sqr present ⇒ oxidative. Annotates the dsrA/dsrB evidence_source
        # (e.g. "ko|dsr_reductive"); the complete_sulfate_reduction vs
        # reverse_dsr_sulfur_oxidation synergies encode the same call downstream.
        coords = parse_prodigal_coords(args.gene_coords)
        by_id = {r["target_id"]: r for r in rows}
        direction = resolve_dsr_direction(by_id, coords)
        if direction:
            for tid in ("dsrA", "dsrB"):
                r = by_id.get(tid)
                if r and r["status"] in ("confirmed", "domain-only", "narrow-no-IPR"):
                    r["evidence_source"] = f"{r['evidence_source'] or 'blast'}|dsr_{direction}"
            print(f"[apply_rules] {args.sample}: dsrAB direction = {direction}"
                  f"{' (synteny-refined)' if coords else ' (presence-based)'}",
                  file=sys.stderr)

        # dsrC / dsrD are dissimilatory-Dsr subunits — an orphan hit with NO dsrA in the
        # genome is the ubiquitous TusE/PF04358 paralog (dsrC) or a stray match (dsrD), not
        # dissimilatory dsr. Demote them to absent. (SP2 fix: cleared 5 dsrC TusE FPs.)
        dsrA_present = (by_id.get("dsrA") or {}).get("status") in (
            "confirmed", "domain-only", "narrow-no-IPR")
        if not dsrA_present:
            for tid in ("dsrC", "dsrD"):
                r = by_id.get(tid)
                if r and r["status"] != "absent":
                    r.update(status="absent", evidence_source="", protein_id="",
                             pfam_hits="", best_pfam_evalue="", blast_acc="",
                             blast_pident="", blast_evalue="")

        # NOTE: an fccA←fccB co-occurrence rescue was tried (SP4b) but REMOVED after the SP6c
        # specialist audit. K17229 ("fccB") is also assigned by KEGG to the Sox-system soxF
        # flavoprotein, so "fccB present" does NOT imply a flavocytochrome-c (FccAB) complex
        # with a cognate FccA cyt subunit (e.g. R. denitrificans: its K17229 hit is soxF in a
        # soxCDEF cluster, and the rescued cyt was an unrelated orphan). FccA (the diheme cyt
        # subunit) cannot be reliably distinguished from Sox / other c-cytochromes without
        # operon synteny, which is unavailable on proteome input — so fccA is left BLAST-gated
        # only (synteny-limited; see config notes + REPORT.md). Genuine FccAB does occur in the
        # panel (e.g. P. denitrificans, operon-verified) but is not reliably callable here.

    elif args.mode == "read":
        rp = load_read_presence(args.read_presence)
        rows = []
        for t in targets:
            res = evaluate_target_from_reads(t, rp)
            if res is None:
                rows.append({
                    "target_id": t["id"], "name": t["name"],
                    "category": t["category"], "complex": t.get("complex", ""),
                    "status": "absent", "evidence_source": "",
                    "protein_id": "", "pfam_hits": "",
                    "best_pfam_evalue": "", "blast_acc": "",
                    "blast_pident": "", "blast_evalue": "",
                    "breadth_pct": rp.get(t["id"], {}).get("breadth_pct", ""),
                    "depth_median": rp.get(t["id"], {}).get("depth_median", ""),
                })
                continue
            _, ev = res
            rows.append({
                "target_id": t["id"], "name": t["name"],
                "category": t["category"], "complex": t.get("complex", ""),
                "status": ev["status"],
                "evidence_source": "read-coverage",
                "protein_id": "", "pfam_hits": "",
                "best_pfam_evalue": "",
                "blast_acc": ev["blast_acc"],
                "blast_pident": ev["blast_pident"], "blast_evalue": "",
                "breadth_pct": ev.get("breadth_pct", ""),
                "depth_median": ev.get("depth_median", ""),
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    n_status = defaultdict(int)
    for r in rows:
        n_status[r["status"]] += 1
    print(f"[apply_rules] {args.sample}: " +
          ", ".join(f"{k}={v}" for k, v in sorted(n_status.items())),
          file=sys.stderr)


if __name__ == "__main__":
    main()
