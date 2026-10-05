#!/usr/bin/env python3
"""
make_html_report.py — one self-contained, interactive HTML report for a run.

Collects, for every sample, the calls (with their evidence), the complex and
process-module completeness, the cycle-map step states and — for nucleotide /
MAG samples — the gene neighbourhoods, and inlines them as JSON into
report_template.html. The page needs no network access and no other file:
it draws the gene grid, the complex / process grid, and a per-genome panel
(cycle map, locus maps, calls table) in the browser, with hover details.
The page font (JetBrains Mono, fonts/) is inlined too, so the report looks the
same offline and its PNG export can use it.
Identical in the nitrogen, sulfur and methane sister pipelines.

Usage:
  python workflow/scripts/make_html_report.py \
      --samples A,B,... --results-dir results \
      --targets config/targets.yaml --out <results>/<tool>_report.html
"""

from __future__ import annotations

import argparse
import base64
import datetime
import json
from pathlib import Path

import yaml

import _cycle_model as model
from _completeness import load_complexes, load_synergies, pretty
from _loci import loci_for_sample
from _viz import (CALLS_TSV, CAT_ORDER, CYCLE_LETTER, CYCLE_NAME, LOCI_TSV,
                  RULED_OUT_LABEL,
                  PATHWAY_COLOR, PATHWAY_COLOR_DARK, PATHWAY_LABEL,
                  PATHWAY_SHORT, PRESENT_CODES, STATUS_CODE, load_calls,
                  order_by_profile)

TEMPLATE = Path(__file__).with_name("report_template.html")
FONT_DIR = Path(__file__).with_name("fonts")
FONTS = (("JetBrainsMono-Regular.woff2", 400), ("JetBrainsMono-Medium.woff2", 500),
         ("JetBrainsMono-Bold.woff2", 700))


def font_faces() -> str:
    """@font-face rules with the report font inlined as base64. A missing file
    is skipped: the page then falls back to the system monospace face."""
    rules = []
    for fname, weight in FONTS:
        path = FONT_DIR / fname
        if not path.exists():
            print(f"[html_report] ! {path} not found — system monospace will be used")
            continue
        b64 = base64.b64encode(path.read_bytes()).decode()
        rules.append('@font-face { font-family: "JetBrains Mono"; font-style: normal; '
                     f'font-weight: {weight}; font-display: block; '
                     f'src: url(data:font/woff2;base64,{b64}) format("woff2"); }}')
    return "\n".join(rules)


def _num(v: str):
    """Keep table numbers as numbers where they are numbers."""
    try:
        return int(v)
    except (TypeError, ValueError):
        return v or ""


def sample_record(name: str, idx: int, results_dir: Path, targets: list[dict],
                  max_gap: int, flank: int) -> dict:
    cdir = results_dir / name / "calls"
    calls = load_calls(cdir / CALLS_TSV)
    codes = {tid: STATUS_CODE.get(r["status"], 0) for tid, r in calls.items()}

    rec_calls = {}
    for tid, r in calls.items():
        if codes[tid] == 0:
            continue
        rec_calls[tid] = [
            codes[tid], r.get("evidence_source", ""), r.get("protein_id", ""),
            r.get("pfam_hits", ""),
            "" if r.get("best_pfam_evalue") == "None" else r.get("best_pfam_evalue", ""),
            r.get("blast_acc", ""), r.get("blast_pident", ""),
            r.get("contig", ""), _num(r.get("start")), _num(r.get("end")),
            r.get("strand", ""),
            len([x for x in (r.get("other_copies") or "").split(";") if x]),
        ]

    cx = {k: [v["state"], v["text"], v["present"], v["missing"]]
          for k, v in load_complexes(cdir / "complex_completeness.tsv").items()}
    syn = {k: [v["state"], v["text"], v["present"], v["missing"], v["violated"]]
           for k, v in load_synergies(cdir / "synergy_completeness.tsv").items()}

    ctx = model.genome_context(calls)
    steps = {sid: [ev["state"],
                   [[[tok["label"], tok["code"]] for tok in line]
                    for line in ev["lines"]],
                   int(ev["reversed"])]
             for sid, ev in model.evaluate(codes, ctx).items()}

    has_coords = (results_dir / name / "prodigal" / f"{name}.faa").exists()
    loci = []
    if has_coords:
        for l in loci_for_sample(results_dir, name, calls, targets, max_gap, flank):
            loci.append({
                "id": l["locus_id"], "pathway": l["pathway"],
                "contig": l["contig"], "clen": l["contig_len"],
                "start": l["start"], "end": l["end"],
                "el": l["edge_left"], "er": l["edge_right"],
                "np": l["n_present"],
                "genes": [[g["id"], g["start"], g["end"], g["strand"],
                           [[c["target"], c["code"], c["category"], int(c["extra"])]
                            for c in g["calls"]], g["partial"]]
                          for g in l["genes"]],
            })
    return {"name": name, "idx": idx,
            "kind": "nucleotide" if has_coords else "protein",
            "note": model.context_note(ctx),
            "calls": rec_calls, "cx": cx, "syn": syn, "steps": steps,
            "loci": loci, "_codes": codes}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="comma-separated names")
    ap.add_argument("--results-dir", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--title", default=f"{CYCLE_NAME.capitalize()}-cycle report")
    ap.add_argument("--max-gap", type=int, default=5)
    ap.add_argument("--flank", type=int, default=2)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    targets = cfg["targets"]
    ordered = [t for cat in CAT_ORDER for t in targets if t["category"] == cat]

    names = [s for s in args.samples.split(",") if s]
    samples = [sample_record(n, i, args.results_dir, targets, args.max_gap, args.flank)
               for i, n in enumerate(names)]

    # Row order "similar gene content" — same rule as the static overview grid.
    vectors = [[1 if s["_codes"].get(t["id"], 0) in PRESENT_CODES else 0
                for t in ordered] for s in samples]
    for rank, i in enumerate(order_by_profile(names, vectors)):
        samples[i]["rank"] = rank
    for s in samples:
        del s["_codes"]

    cats = [c for c in CAT_ORDER if any(t["category"] == c for t in targets)]
    data = {
        "generated": datetime.date.today().isoformat(),
        "labels": {"cycle": CYCLE_NAME, "letter": CYCLE_LETTER,
                   "loci_tsv": LOCI_TSV, "ruled_out": RULED_OUT_LABEL},
        "pathways": [{"id": c, "label": PATHWAY_LABEL[c], "short": PATHWAY_SHORT[c]}
                     for c in cats],
        "targets": [{"id": t["id"], "label": model.DISPLAY.get(t["id"], t["id"]),
                     "name": t["name"], "cat": t["category"]} for t in ordered],
        "complexes": [{"id": cid, "label": pretty(cid),
                       "desc": (body.get("application", "").split(";")[0]).strip()}
                      for cid, body in cfg.get("complexes", {}).items()],
        "synergies": [{"id": s["name"], "label": pretty(s["name"]),
                       "desc": s.get("benefit", "").split(". ")[0].strip().rstrip("."),
                       "group": "flag" if s.get("forbids") else "process"}
                      for s in cfg.get("synergies", [])],
        "cycle": model.layout(),
        "samples": samples,
    }

    css_light = " ".join(f"--p{i}: {PATHWAY_COLOR[c]};" for i, c in enumerate(cats))
    css_dark = " ".join(f"--p{i}: {PATHWAY_COLOR_DARK[c]};" for i, c in enumerate(cats))
    # `</` must not appear inside the inline <script> block.
    blob = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")

    html = TEMPLATE.read_text()
    for key, val in (("/*__FONT_FACES__*/", font_faces()),
                     ("/*__PATHWAY_CSS_LIGHT__*/", css_light),
                     ("/*__PATHWAY_CSS_DARK__*/", css_dark),
                     ("__TITLE__", args.title),
                     ("__FILE__", args.out.name),
                     ("__DATA_JSON__", blob)):
        if key not in html:
            raise SystemExit(f"report template is missing the {key} placeholder")
        html = html.replace(key, val)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    print(f"[html_report] wrote {args.out}  ({len(samples)} samples, "
          f"{len(html) / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
