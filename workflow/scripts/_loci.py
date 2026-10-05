"""Gene-neighbourhood ("locus") assembly for nucleotide / MAG samples.

Reads the Prodigal proteome (coordinates live in the FASTA headers) and the
per-sample calls, clusters the called genes that sit close together on a
contig, and returns each cluster with its flanking genes. Used by
make_locus_maps.py (static figure + the loci table) and make_html_report.py.

Only genes in the target set are annotated; every other gene in a locus window
is returned with no call (drawn blank).
"""

from __future__ import annotations

import re
from pathlib import Path

from _viz import CAT_ORDER, STATUS_CODE

# Rank used to pick the gene label when one protein carries several calls
# (e.g. amoA confirmed + amoA_gamma disqualified on the same ORF).
_RANK = {2: 0, 1: 1, -1: 2, 0: 3}


def parse_prodigal_genes(faa: Path) -> dict[str, list[dict]]:
    """{contig: [gene, ...]} in genomic order. gene = {id, contig, idx, start,
    end, strand, partial}. Empty dict when the FASTA has no Prodigal headers."""
    contigs: dict[str, list[dict]] = {}
    if not faa or not Path(faa).exists():
        return contigs
    with open(faa) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = [p.strip() for p in line[1:].split("#")]
            if len(parts) < 4:
                continue
            pid = parts[0].split()[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            contig = pid.rsplit("_", 1)[0]
            m = re.search(r"partial=(\d\d)", parts[4] if len(parts) > 4 else "")
            contigs.setdefault(contig, []).append({
                "id": pid, "contig": contig, "start": min(start, end),
                "end": max(start, end),
                "strand": "-" if parts[3] == "-1" else "+",
                "partial": m.group(1) if m else "00",
            })
    for genes in contigs.values():
        genes.sort(key=lambda g: g["start"])
        for i, g in enumerate(genes):
            g["idx"] = i
    return contigs


def parse_contig_lengths(gff: Path) -> dict[str, int]:
    """{contig: length} from the `# Sequence Data:` lines of a Prodigal GFF."""
    out: dict[str, int] = {}
    if not gff or not Path(gff).exists():
        return out
    with open(gff) as fh:
        for line in fh:
            if not line.startswith("# Sequence Data:"):
                continue
            m_len = re.search(r"seqlen=(\d+)", line)
            m_hdr = re.search(r'seqhdr="([^"]*)"', line)
            if m_len and m_hdr and m_hdr.group(1).split():
                out[m_hdr.group(1).split()[0]] = int(m_len.group(1))
    return out


def build_loci(calls: dict[str, dict], contigs: dict[str, list[dict]],
               contig_len: dict[str, int], targets: list[dict],
               max_gap: int = 5, flank: int = 2) -> list[dict]:
    """Cluster called genes into loci.

    Two called genes on one contig join the same locus when at most `max_gap`
    other genes lie between them. Each locus is returned with `flank` extra
    genes on either side. A locus is filed under the pathway that contributes
    most of its present (confirmed / domain-only) calls.
    """
    cat_of = {t["id"]: t["category"] for t in targets}
    by_id = {g["id"]: g for genes in contigs.values() for g in genes}

    # protein id → list of (target_id, code, is_extra_copy), best call first.
    # The calls table reports one protein per target; proteins listed under
    # `other_copies` reach the same status and are drawn too (flagged), so a
    # second operon copy does not show up as a row of blank arrows.
    hits: dict[str, list[tuple[str, int, bool]]] = {}
    for tid, r in calls.items():
        code = STATUS_CODE.get(r["status"], 0)
        if code == 0:
            continue
        if r.get("protein_id") in by_id:
            hits.setdefault(r["protein_id"], []).append((tid, code, False))
        if code in (1, 2):
            for pid in (r.get("other_copies") or "").split(";"):
                if pid in by_id:
                    hits.setdefault(pid, []).append((tid, code, True))
    for lst in hits.values():
        lst.sort(key=lambda x: (_RANK[x[1]], x[2], x[0]))

    loci: list[dict] = []
    for contig, genes in contigs.items():
        idxs = sorted(by_id[p]["idx"] for p in hits if by_id[p]["contig"] == contig)
        if not idxs:
            continue
        groups = [[idxs[0]]]
        for i in idxs[1:]:
            if i - groups[-1][-1] - 1 <= max_gap:
                groups[-1].append(i)
            else:
                groups.append([i])
        for grp in groups:
            lo = max(0, grp[0] - flank)
            hi = min(len(genes) - 1, grp[-1] + flank)
            window = []
            votes: dict[str, int] = {}
            n_present = 0
            for g in genes[lo:hi + 1]:
                g_hits = hits.get(g["id"], [])
                item = dict(g)
                item["calls"] = [{"target": t, "code": c, "extra": x,
                                  "category": cat_of.get(t, "")}
                                 for t, c, x in g_hits]
                window.append(item)
                if g_hits and g_hits[0][1] in (1, 2):
                    n_present += 1
                    cat = cat_of.get(g_hits[0][0], "")
                    votes[cat] = votes.get(cat, 0) + 1
            if not votes:      # only disqualified hits: file under the first one
                first = next(i for i in window if i["calls"])
                votes[first["calls"][0]["category"]] = 0
            pathway = min(votes, key=lambda c: (-votes[c],
                                                CAT_ORDER.index(c) if c in CAT_ORDER else 99))
            length = contig_len.get(contig, 0)
            loci.append({
                "contig": contig,
                "contig_len": length,
                "start": window[0]["start"],
                "end": window[-1]["end"],
                "genes": window,
                "n_present": n_present,
                "pathway": pathway,
                # The window touches the first / last gene of the contig: the
                # operon may continue beyond the assembly here.
                "edge_left": lo == 0,
                "edge_right": hi == len(genes) - 1,
            })
    loci.sort(key=lambda l: (CAT_ORDER.index(l["pathway"]) if l["pathway"] in CAT_ORDER else 99,
                             -l["n_present"], l["contig"], l["start"]))
    for k, l in enumerate(loci, 1):
        l["locus_id"] = f"L{k:03d}"
    return loci


def loci_for_sample(results_dir: Path, sample: str, calls: dict[str, dict],
                    targets: list[dict], max_gap: int = 5, flank: int = 2) -> list[dict]:
    """Loci for one sample, or [] when it has no Prodigal output (protein input)."""
    pdir = Path(results_dir) / sample / "prodigal"
    contigs = parse_prodigal_genes(pdir / f"{sample}.faa")
    if not contigs:
        return []
    lengths = parse_contig_lengths(pdir / f"{sample}.gff")
    return build_loci(calls, contigs, lengths, targets, max_gap, flank)
