"""The sulfur-cycle diagram: compounds, reaction steps, and how a genome's
calls light the steps up.

Single source of truth for the static cycle maps (make_cycle_map.py) and
the interactive one in report.html (make_html_report.py ships the geometry
computed here as JSON, so the browser only draws polylines).

A step is `complete` when every gene of at least one route is present
(confirmed or domain-only — the same rule as complex completeness), `partial`
when some route gene is present but no route is whole, otherwise `absent`.
This is a drawing-level summary; the authoritative per-process calls remain
complex_completeness.tsv and synergy_completeness.tsv.

The sat / aprAB / dsrAB chain is reversible: it is drawn SO₄²⁻ → H₂S unless
apply_rules.py called the genome's dsrAB `oxidative` (reverse Dsr), in which
case the three arrows point the other way.
"""

from __future__ import annotations

import math
import re

# Canvas (data units, y up).
XLIM = (-3.0, 181.0)
YLIM = (-3.0, 106.0)

# id: (label, x, y, width, height)
NODES = {
    "SO4":  ("SO₄²⁻",       8.0, 66.0, 13.0, 7.0),
    "APS":  ("APS",         36.0, 90.0, 10.0, 7.0),
    "PAPS": ("PAPS",        86.0, 97.0, 11.0, 7.0),
    "SO3":  ("SO₃²⁻",       78.0, 66.0, 13.0, 7.0),
    "H2S":  ("H₂S",        120.0, 66.0, 10.0, 7.0),
    "S0":   ("S⁰",          99.0, 40.0,  8.0, 7.0),
    "SULF": ("sulfonates",  72.0, 41.0, 19.0, 7.0),
    "S2O3": ("S₂O₃²⁻",      40.0, 18.0, 15.0, 7.0),
    "S4O6": ("S₄O₆²⁻",       8.0, 18.0, 15.0, 7.0),
    "CYS":  ("cysteine",   152.0, 66.0, 16.0, 7.0),
    "MESH": ("MeSH",       144.0, 91.0, 12.0, 7.0),
    "DMSP": ("DMSP",       170.0, 99.0, 12.0, 7.0),
    "DMS":  ("DMS",        170.0, 77.0, 10.0, 7.0),
}

# Short display names (none needed for the sulfur targets).
DISPLAY: dict[str, str] = {}

_D = "dissimilatory_sulfate_reduction"
_O = "sulfur_oxidation"
_A = "assimilatory_sulfate_reduction"
_T = "thiosulfate_polysulfide"
_M = "organic_sulfur_dmsp"
_U = "sulfonate_taurine"

# arrows: (from, to, bow, (dx0, dy0), (dx1, dy1)) — bow is the offset of the
#   curve's control point to the LEFT of the direction of travel; the two
#   offsets shift the start / end anchor off the node centre (parallel arrows).
# routes: alternative gene sets, any one of which performs the step.
# lines:  how the genes are written next to the arrow (may include accessory
#   genes that no route requires, e.g. qmoABC, dsrC/D/M/K, tauABC).
# label:  (x, y, ha) of the first label line; further lines stack downwards.
# reversible: drawn the other way round in a reverse-Dsr (oxidative) genome.
STEPS = [
    {"id": "sat", "name": "Sulfate activation / APS → sulfate (ATP sulfurylase)",
     "pathway": _D, "reversible": True,
     "arrows": [("SO4", "APS", 0.0, (-1.6, 1.4), (-1.6, 1.4))],
     "routes": [["sat"]], "lines": [["sat"]],
     "label": (17.0, 82.0, "right")},
    {"id": "cysDN", "name": "Sulfate activation (assimilatory)", "pathway": _A,
     "arrows": [("SO4", "APS", 0.0, (1.6, -1.4), (1.6, -1.4))],
     "routes": [["cysN", "cysD"], ["cysNC"]],
     "lines": [["cysN", "cysD"], ["cysNC"]],
     "label": (26.0, 73.5, "left")},
    {"id": "apr", "name": "APS ⇄ sulfite (APS reductase)", "pathway": _D,
     "reversible": True,
     "arrows": [("APS", "SO3", 0.0, (0, 0), (0, 0))],
     "routes": [["aprA", "aprB"]],
     "lines": [["aprA", "aprB"], ["qmoA", "qmoB", "qmoC"]],
     "label": (61.0, 86.0, "left")},
    {"id": "cysC", "name": "APS → PAPS (APS kinase)", "pathway": _A,
     "arrows": [("APS", "PAPS", -3.0, (0, 1.0), (0, 0))],
     "routes": [["cysC"], ["cysNC"]], "lines": [["cysC"]],
     "label": (60.0, 100.5, "center")},
    {"id": "cysH", "name": "PAPS → sulfite (PAPS / APS reductase)", "pathway": _A,
     "arrows": [("PAPS", "SO3", 0.0, (0, 0), (2.0, 0))],
     "routes": [["cysH"]], "lines": [["cysH"]],
     "label": (86.5, 84.0, "left")},
    {"id": "dsr", "name": "Sulfite ⇄ sulfide (dissimilatory sulfite reductase)",
     "pathway": _D, "reversible": True,
     "arrows": [("SO3", "H2S", 0.0, (0, 1.6), (0, 1.6))],
     "routes": [["dsrA", "dsrB"]],
     "lines": [["dsrA", "dsrB"], ["dsrC", "dsrD"], ["dsrM", "dsrK"]],
     "label": (99.0, 80.5, "center")},
    {"id": "sir", "name": "Sulfite → sulfide (assimilatory sulfite reductase)",
     "pathway": _A,
     "arrows": [("SO3", "H2S", 0.0, (0, -1.6), (0, -1.6))],
     "routes": [["cysJ", "cysI"], ["sir"]],
     "lines": [["cysJ", "cysI"], ["sir"]],
     "label": (99.0, 59.5, "center")},
    {"id": "cysK", "name": "Sulfide → cysteine (cysteine synthase)", "pathway": _A,
     "arrows": [("H2S", "CYS", 0.0, (0, 0), (0, 0))],
     "routes": [["cysK"], ["cysM"]], "lines": [["cysK", "cysM"]],
     "label": (135.0, 61.0, "center")},
    {"id": "sqr", "name": "Sulfide → S⁰ (SQR / flavocytochrome c)", "pathway": _O,
     "arrows": [("H2S", "S0", 0.0, (-1.6, 1.2), (-1.6, 1.2))],
     "routes": [["sqr"], ["fccA", "fccB"]],
     "lines": [["sqr", "fccA", "fccB"]],
     "label": (99.0, 32.5, "center")},
    {"id": "sre", "name": "S⁰ → sulfide (sulfur reductase)", "pathway": _T,
     "arrows": [("S0", "H2S", 0.0, (1.6, -1.2), (1.6, -1.2))],
     "routes": [["sreA"]], "lines": [["sreA"]],
     "label": (112.5, 50.0, "left")},
    {"id": "sdo", "name": "S⁰ → sulfite (sulfur dioxygenase / SOR)", "pathway": _O,
     "arrows": [("S0", "SO3", 0.0, (0, 0), (3.5, 0))],
     "routes": [["sdo"], ["sor"]], "lines": [["sdo", "sor"]],
     "label": (86.0, 50.5, "right")},
    {"id": "soe", "name": "Sulfite → sulfate (sulfite dehydrogenase)", "pathway": _O,
     "arrows": [("SO3", "SO4", 0.0, (0, 0), (0, 0))],
     "routes": [["soeA", "soeB", "soeC"], ["sorA"]],
     "lines": [["soeA", "soeB", "soeC", "sorA"]],
     "label": (43.0, 61.5, "center")},
    {"id": "sox", "name": "Thiosulfate → sulfate (Sox system)", "pathway": _O,
     "arrows": [("S2O3", "SO4", 2.0, (-2.0, 0), (0, 0))],
     "routes": [["soxA", "soxB", "soxX", "soxY", "soxZ"]],
     "lines": [["soxA", "soxB", "soxX"], ["soxY", "soxZ"], ["soxC", "soxD"]],
     "label": (24.0, 51.0, "left")},
    {"id": "tsd", "name": "Thiosulfate → tetrathionate (TsdA)", "pathway": _O,
     "arrows": [("S2O3", "S4O6", 0.0, (0, 2.3), (0, 2.3))],
     "routes": [["tsdA"]], "lines": [["tsdA"]],
     "label": (24.0, 25.5, "center")},
    {"id": "dox", "name": "Thiosulfate → tetrathionate (TQO, DoxDA)", "pathway": _T,
     "arrows": [("S2O3", "S4O6", 0.0, (0, 0), (0, 0))],
     "routes": [["doxA", "doxD"]], "lines": [["doxA", "doxD"]],
     "label": (24.0, 11.0, "center")},
    {"id": "ttr", "name": "Tetrathionate → thiosulfate (Ttr / Otr / Tth)",
     "pathway": _T,
     "arrows": [("S4O6", "S2O3", 0.0, (0, -2.3), (0, -2.3))],
     "routes": [["ttrA", "ttrB", "ttrC"], ["otr"], ["tth"]],
     "lines": [["ttrA", "ttrB", "ttrC"], ["otr", "tth"]],
     "label": (24.0, 7.2, "center")},
    {"id": "phs", "name": "Thiosulfate → sulfide + sulfite (thiosulfate reductase)",
     "pathway": _T,
     "arrows": [("S2O3", "H2S", -60.0, (0, -1.0), (4.0, -1.0))],
     "routes": [["phsA", "phsB", "phsC"]],
     "lines": [["phsA", "phsB", "phsC"]],
     "label": (88.0, 6.0, "center")},
    {"id": "tau", "name": "Sulfonate desulfonation → sulfite", "pathway": _U,
     "arrows": [("SULF", "SO3", 0.0, (1.0, 0), (0.5, 0))],
     "routes": [["tauD"], ["ssuD"]],
     "lines": [["tauA", "tauB", "tauC", "tauD"], ["ssuD", "ssuE"]],
     "label": (72.0, 33.5, "center")},
    {"id": "sse", "name": "Thiosulfate → sulfite (rhodanese, thiosulfate sulfurtransferase)",
     "pathway": _T,
     "arrows": [("S2O3", "SO3", 5.0, (3.0, 0), (-4.0, 0))],
     "routes": [["sseA"]], "lines": [["sseA"]],
     "label": (47.5, 34.5, "right")},
    {"id": "mto", "name": "Methanethiol → sulfide (methanethiol oxidase)",
     "pathway": _M,
     "arrows": [("MESH", "H2S", 0.0, (0, 0), (1.0, 0))],
     "routes": [["mtoX"]], "lines": [["mtoX"]],
     "label": (136.0, 77.0, "left")},
    {"id": "dmd", "name": "DMSP demethylation → methanethiol", "pathway": _M,
     "arrows": [("DMSP", "MESH", 0.0, (0, 0), (0, 0))],
     "routes": [["dmdA"]], "lines": [["dmdA"]],
     "label": (156.0, 99.5, "center")},
    {"id": "ddd", "name": "DMSP cleavage → DMS (DMSP lyase)", "pathway": _M,
     "arrows": [("DMSP", "DMS", 0.0, (0, 0), (0, 0))],
     "routes": [["dddP"]], "lines": [["dddP"]],
     "label": (172.5, 88.0, "left")},
    {"id": "mdd", "name": "Methanethiol → DMS (MddA)", "pathway": _M,
     "arrows": [("MESH", "DMS", 0.0, (0, 0), (0, 0))],
     "routes": [["mddA"]], "lines": [["mddA"]],
     "label": (155.0, 80.0, "right")},
]

# Figure sizing used by make_cycle_map.py (inches / points).
FIG = {"single_w": 11.4, "grid_w": 8.4, "single_fs": 9.0, "grid_fs": 6.6,
       "grid_cols": 2}

LINE_STEP = 3.6          # distance between stacked label lines (data units)
HEAD_LEN, HEAD_W = 2.4, 1.25
NODE_PAD = 1.3           # gap between an arrow end and the node box


# ───────────────────────────── evaluation ────────────────────────────────────

def _slot_genes(slot) -> list[str]:
    return [slot] if isinstance(slot, str) else list(slot)


def _token(slot, codes: dict[str, int]) -> dict:
    """The gene to print for a slot: the best-supported alternative (first one
    when none is present)."""
    genes = _slot_genes(slot)
    rank = {2: 0, 1: 1, -1: 2, 0: 3}
    best = min(genes, key=lambda g: (rank[codes.get(g, 0)], genes.index(g)))
    return {"gene": best, "label": DISPLAY.get(best, best),
            "code": codes.get(best, 0)}


def genome_context(calls: dict[str, dict]) -> dict:
    """dsrAB runs in either direction and sequence alone cannot tell which;
    apply_rules.py calls it from the genomic companions (dsrD / qmo ⇒ reductive,
    Sox / sqr ⇒ oxidative) and tags it on the dsrA / dsrB evidence_source, e.g.
    `ko|dsr_oxidative`. Read that tag back."""
    for tid in ("dsrA", "dsrB"):
        m = re.search(r"dsr_(reductive|oxidative|ambiguous)",
                      (calls.get(tid) or {}).get("evidence_source", ""))
        if m:
            return {"dsr_direction": m.group(1)}
    return {}


def context_note(ctx: dict) -> str:
    """One-line description of the context, for subtitles ('' if none)."""
    d = (ctx or {}).get("dsr_direction")
    if d == "oxidative":
        return "dsrAB direction: oxidative (reverse Dsr)"
    if d == "reductive":
        return "dsrAB direction: reductive"
    if d == "ambiguous":
        return "dsrAB direction: ambiguous (drawn reductive)"
    return ""


def evaluate(codes: dict[str, int], ctx: dict | None = None) -> dict[str, dict]:
    """{step_id: {state, reversed, lines:[[token,…],…]}} for one genome."""
    out: dict[str, dict] = {}
    for st in STEPS:
        def present(slot) -> bool:
            return any(codes.get(g, 0) in (1, 2) for g in _slot_genes(slot))
        route_full = any(all(present(s) for s in r) for r in st["routes"])
        any_gene = any(present(s) for r in st["routes"] for s in r)
        state = "complete" if route_full else "partial" if any_gene else "absent"
        out[st["id"]] = {
            "state": state,
            "reversed": bool(st.get("reversible")
                             and (ctx or {}).get("dsr_direction") == "oxidative"),
            "lines": [[_token(s, codes) for s in line] for line in st["lines"]],
        }
    return out


# ───────────────────────────── geometry ──────────────────────────────────────

def _inside(pt, node, pad) -> bool:
    _, x, y, w, h = NODES[node]
    return abs(pt[0] - x) <= w / 2 + pad and abs(pt[1] - y) <= h / 2 + pad


def arrow_geometry(frm, to, bow, off0, off1, n=48) -> dict:
    """Polyline (trimmed at both node boxes) + arrowhead triangle."""
    p0 = (NODES[frm][1] + off0[0], NODES[frm][2] + off0[1])
    p2 = (NODES[to][1] + off1[0], NODES[to][2] + off1[1])
    dx, dy = p2[0] - p0[0], p2[1] - p0[1]
    d = math.hypot(dx, dy) or 1.0
    c = ((p0[0] + p2[0]) / 2 - bow * dy / d, (p0[1] + p2[1]) / 2 + bow * dx / d)
    pts = []
    for i in range(n + 1):
        t = i / n
        a, b, e = (1 - t) ** 2, 2 * (1 - t) * t, t ** 2
        pts.append((a * p0[0] + b * c[0] + e * p2[0],
                    a * p0[1] + b * c[1] + e * p2[1]))
    pts = [p for p in pts
           if not _inside(p, frm, NODE_PAD) and not _inside(p, to, NODE_PAD)]
    if len(pts) < 2:
        return {"line": [], "head": []}
    tip = pts[-1]
    # Direction at the tip, taken a few samples back for stability.
    ref = pts[max(0, len(pts) - 4)]
    ux, uy = tip[0] - ref[0], tip[1] - ref[1]
    u = math.hypot(ux, uy) or 1.0
    ux, uy = ux / u, uy / u
    base = (tip[0] - ux * HEAD_LEN, tip[1] - uy * HEAD_LEN)
    head = [tip,
            (base[0] - uy * HEAD_W, base[1] + ux * HEAD_W),
            (base[0] + uy * HEAD_W, base[1] - ux * HEAD_W)]
    # End the shaft at the arrowhead base so a thick line never pokes through.
    shaft = [p for p in pts
             if (p[0] - tip[0]) * ux + (p[1] - tip[1]) * uy <= -HEAD_LEN * 0.8]
    shaft.append(base)
    return {"line": shaft, "head": head}


def layout() -> dict:
    """Everything a renderer needs, JSON-serialisable."""
    r = lambda p: [round(p[0], 2), round(p[1], 2)]          # noqa: E731
    steps = []
    for st in STEPS:
        def geom(specs):
            out = []
            for a in specs:
                g = arrow_geometry(*a)
                out.append({"line": [r(p) for p in g["line"]],
                            "head": [r(p) for p in g["head"]]})
            return out
        # A reversible step also carries the same arrows drawn the other way.
        rev = [(to, frm, -bow, o1, o0) for frm, to, bow, o0, o1 in st["arrows"]] \
            if st.get("reversible") else []
        steps.append({"id": st["id"], "name": st["name"],
                      "pathway": st["pathway"], "arrows": geom(st["arrows"]),
                      "arrows_rev": geom(rev), "label": list(st["label"])})
    return {
        "xlim": list(XLIM), "ylim": list(YLIM), "line_step": LINE_STEP,
        "wide": (XLIM[1] - XLIM[0]) > 130,
        "nodes": [{"id": k, "label": v[0], "x": v[1], "y": v[2],
                   "w": v[3], "h": v[4]} for k, v in NODES.items()],
        "steps": steps,
    }
