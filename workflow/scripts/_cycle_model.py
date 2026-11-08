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

The diagram is drawn as a metro map: each pathway is a line in its colour,
each compound a station, each gene a stop on the line of the step it performs.
Sulfate reduction is the trunk from SO₄²⁻ (left) over APS to SO₃²⁻, H₂S and
cysteine (right), with the PAPS detour above APS → SO₃²⁻ and the dissimilatory
and assimilatory lines as parallel strands; sulfite oxidation (Soe) returns
along the trunk; S⁰, sulfonates, thiosulfate and tetrathionate sit below it,
the DMSP / methanethiol / DMS triangle above H₂S.

The sat / aprAB / dsrAB chain is reversible: it is travelled SO₄²⁻ → H₂S unless
apply_rules.py called the genome's dsrAB `oxidative` (reverse Dsr), in which
case the three lines run the other way.
"""

from __future__ import annotations

import math
import re

# Canvas (data units, y up).
XLIM = (0.0, 276.0)
YLIM = (0.0, 126.0)

# Stations.  id: (label, x, y, (label dx, label dy, ha, va))
# Lines meet at the station centre; the label offset keeps the name clear of
# the lines that arrive there.
NODES = {
    "SO4":  ("SO₄²⁻",       34.0,  70.0, (-6.5,  0.0, "right",  "center")),
    "APS":  ("APS",         82.0, 100.0, (0.0,  -6.0, "center", "top")),
    "PAPS": ("PAPS",       112.0, 118.0, (6.5,   0.0, "left",   "center")),
    "SO3":  ("SO₃²⁻",      142.0,  70.0, (5.0,   5.0, "left",   "bottom")),
    "H2S":  ("H₂S",        202.0,  70.0, (5.0,   5.0, "left",   "bottom")),
    "S0":   ("S⁰",         238.0,  34.0, (6.5,   0.0, "left",   "center")),
    "SULF": ("sulfonates", 142.0,  14.0, (-6.5,  0.0, "right",  "center")),
    "S2O3": ("S₂O₃²⁻",      82.0,  22.0, (-1.5, -9.5, "center", "top")),
    "S4O6": ("S₄O₆²⁻",      26.0,  22.0, (-6.5,  0.0, "right",  "center")),
    "CYS":  ("cysteine",   250.0,  70.0, (6.5,   0.0, "left",   "center")),
    "MESH": ("MeSH",       202.0, 106.0, (-6.5,  0.0, "right",  "center")),
    "DMSP": ("DMSP",       250.0, 106.0, (6.5,   0.0, "left",   "center")),
    "DMS":  ("DMS",        226.0,  88.0, (0.0,   6.0, "center", "bottom")),
}

# Short display names (none needed for the sulfur targets).
DISPLAY: dict[str, str] = {}

_D = "dissimilatory_sulfate_reduction"
_O = "sulfur_oxidation"
_A = "assimilatory_sulfate_reduction"
_T = "thiosulfate_polysulfide"
_M = "organic_sulfur_dmsp"
_U = "sulfonate_taurine"

# path:   waypoints of the line from the first station to the last, octilinear
#   (horizontal / vertical / 45°); the stations' centres are the two ends.
# offset: shift of the whole path to the left of the direction of travel, for
#   parallel strands between the same two stations (0 = on the centre line).
# side:   where the stop labels go — +1 left of the direction of travel, -1
#   right of it, 0 alternating (dense steps); a list gives one value per
#   straight stretch of the path.
# feeders: extra lines of the same step that bring a second substrate to the
#   product station (waypoint lists, drawn in the step's colour and state, with
#   a chevron and no stops; never reversed).
# flat:   write the stop labels horizontally even on a horizontal stretch (a
#   lone stop between parallel strands, where 45° text would cross them).
# stop_shift: moves the group of stops along the path (units, towards the end).
# stop_spacing: this step's own cap on the distance between stops (default
#   STOP_SPACING_MAX).
# routes: alternative gene sets, any one of which performs the step.
# lines:  the genes drawn as stops, in order (may include accessory genes that
#   no route requires, e.g. qmoABC, dsrC/D/M/K, tauABC).
# reversible: travelled the other way round in a reverse-Dsr (oxidative) genome
#   (the chevron moves to the other end).
STEPS = [
    {"id": "sat", "name": "Sulfate activation / APS → sulfate (ATP sulfurylase)",
     "pathway": _D, "reversible": True,
     "path": ["SO4", (64.0, 100.0), "APS"], "offset": 1.6, "side": 1,
     "routes": [["sat"]], "lines": [["sat"]]},
    {"id": "cysDN", "name": "Sulfate activation (assimilatory)", "pathway": _A,
     "path": ["SO4", (64.0, 100.0), "APS"], "offset": -1.6, "side": -1,
     "routes": [["cysN", "cysD"], ["cysNC"]],
     "lines": [["cysN", "cysD"], ["cysNC"]]},
    {"id": "apr", "name": "APS ⇄ sulfite (APS reductase)", "pathway": _D,
     "reversible": True,
     "path": ["APS", (112.0, 100.0), "SO3"], "side": -1,
     "routes": [["aprA", "aprB"]],
     "lines": [["aprA", "aprB"], ["qmoA", "qmoB", "qmoC"]]},
    {"id": "cysC", "name": "APS → PAPS (APS kinase)", "pathway": _A,
     "path": ["APS", (100.0, 118.0), "PAPS"], "side": 1,
     "routes": [["cysC"], ["cysNC"]], "lines": [["cysC"]]},
    {"id": "cysH", "name": "PAPS → sulfite (PAPS / APS reductase)", "pathway": _A,
     "path": ["PAPS", (142.0, 88.0), "SO3"], "side": 1,
     "routes": [["cysH"]], "lines": [["cysH"]]},
    {"id": "dsr", "name": "Sulfite ⇄ sulfide (dissimilatory sulfite reductase)",
     "pathway": _D, "reversible": True,
     "path": ["SO3", "H2S"], "offset": 1.6, "side": 1, "stop_shift": 9.0, "stop_spacing": 7.0,
     "routes": [["dsrA", "dsrB"]],
     "lines": [["dsrA", "dsrB"], ["dsrC", "dsrD"], ["dsrM", "dsrK"]]},
    {"id": "sir", "name": "Sulfite → sulfide (assimilatory sulfite reductase)",
     "pathway": _A,
     "path": ["SO3", "H2S"], "offset": -1.6, "side": -1, "stop_shift": 6.0,
     "routes": [["cysJ", "cysI"], ["sir"]],
     "lines": [["cysJ", "cysI"], ["sir"]]},
    {"id": "cysK", "name": "Sulfide → cysteine (cysteine synthase)", "pathway": _A,
     "path": ["H2S", "CYS"], "side": 1,
     "routes": [["cysK"], ["cysM"]], "lines": [["cysK", "cysM"]]},
    {"id": "sqr", "name": "Sulfide → S⁰ (SQR / flavocytochrome c)", "pathway": _O,
     "path": ["H2S", "S0"], "offset": 1.6, "side": 1,
     "routes": [["sqr"], ["fccA", "fccB"]],
     "lines": [["sqr", "fccA", "fccB"]]},
    {"id": "sre", "name": "S⁰ → sulfide (sulfur reductase)", "pathway": _T,
     "path": ["S0", "H2S"], "offset": 1.6, "side": 1,
     "routes": [["sreA"]], "lines": [["sreA"]]},
    {"id": "sdo", "name": "S⁰ → sulfite (sulfur dioxygenase / SOR)", "pathway": _O,
     "path": ["S0", (178.0, 34.0), "SO3"], "side": 1,
     "routes": [["sdo"], ["sor"]], "lines": [["sdo", "sor"]]},
    {"id": "soe", "name": "Sulfite → sulfate (sulfite dehydrogenase)", "pathway": _O,
     "path": ["SO3", "SO4"], "side": -1,
     "routes": [["soeA", "soeB", "soeC"], ["sorA"]],
     "lines": [["soeA", "soeB", "soeC", "sorA"]]},
    {"id": "sox", "name": "Thiosulfate → sulfate (Sox system)", "pathway": _O,
     "path": ["S2O3", "SO4"], "side": -1, "stop_spacing": 7.5,
     "routes": [["soxA", "soxB", "soxX", "soxY", "soxZ"]],
     "lines": [["soxA", "soxB", "soxX"], ["soxY", "soxZ"], ["soxC", "soxD"]]},
    {"id": "tsd", "name": "Thiosulfate → tetrathionate (TsdA)", "pathway": _O,
     "path": ["S2O3", "S4O6"], "side": -1, "flat": True,
     "routes": [["tsdA"]], "lines": [["tsdA"]]},
    {"id": "dox", "name": "Thiosulfate → tetrathionate (TQO, DoxDA)", "pathway": _T,
     # fans in to a close parallel of the TsdA strand at the S₂O₃²⁻ disc: a 45°
     # approach from the upper left would lie on the Sox line.
     "path": [(82.0, 25.2), (73.8, 25.2), (70.0, 29.0), (33.0, 29.0), "S4O6"], "side": -1,
     "stop_shift": 10.0,
     "routes": [["doxA", "doxD"]], "lines": [["doxA", "doxD"]]},
    {"id": "ttr", "name": "Tetrathionate → thiosulfate (Ttr / Otr / Tth)",
     "pathway": _T,
     "path": ["S4O6", (33.0, 15.0), (75.0, 15.0), "S2O3"], "side": -1, "stop_shift": -4.0,
     "stop_spacing": 6.0,
     "routes": [["ttrA", "ttrB", "ttrC"], ["otr"], ["tth"]],
     "lines": [["ttrA", "ttrB", "ttrC"], ["otr", "tth"]]},
    {"id": "phs", "name": "Thiosulfate → sulfide + sulfite (thiosulfate reductase)",
     "pathway": _T,
     "path": ["S2O3", (98.0, 6.0), (190.0, 6.0), (202.0, 18.0), "H2S"], "side": 1,
     "stop_shift": 8.0,
     "routes": [["phsA", "phsB", "phsC"]],
     "lines": [["phsA", "phsB", "phsC"]]},
    {"id": "tau", "name": "Sulfonate desulfonation → sulfite", "pathway": _U,
     "path": ["SULF", "SO3"], "side": -1, "stop_shift": -4.0, "stop_spacing": 6.5,
     "routes": [["tauD"], ["ssuD"]],
     "lines": [["tauA", "tauB", "tauC", "tauD"], ["ssuD", "ssuE"]]},
    {"id": "sse", "name": "Thiosulfate → sulfite (rhodanese, thiosulfate sulfurtransferase)",
     "pathway": _T,
     "path": ["S2O3", (94.0, 22.0), "SO3"], "side": 1,
     "routes": [["sseA"]], "lines": [["sseA"]]},
    {"id": "mto", "name": "Methanethiol → sulfide (methanethiol oxidase)",
     "pathway": _M,
     "path": ["MESH", "H2S"], "side": 1,
     "routes": [["mtoX"]], "lines": [["mtoX"]]},
    {"id": "dmd", "name": "DMSP demethylation → methanethiol", "pathway": _M,
     "path": ["DMSP", "MESH"], "side": -1,
     "routes": [["dmdA"]], "lines": [["dmdA"]]},
    {"id": "ddd", "name": "DMSP cleavage → DMS (DMSP lyase)", "pathway": _M,
     "path": ["DMSP", (232.0, 88.0), "DMS"], "side": 1,
     "routes": [["dddP"]], "lines": [["dddP"]]},
    {"id": "mdd", "name": "Methanethiol → DMS (MddA)", "pathway": _M,
     "path": ["MESH", (220.0, 88.0), "DMS"], "side": 1,
     "routes": [["mddA"]], "lines": [["mddA"]]},
]

# Figure sizing used by make_cycle_map.py (inches / points).
FIG = {"single_w": 13.2, "grid_w": 9.7, "single_fs": 9.0, "grid_fs": 6.6,
       "grid_cols": 2}

# Metro geometry (data units).
LINE_W = 1.5             # stroke width of a line
STATION_R = 2.5          # station disc radius (ring drawn on it)
STOP_R = 1.05            # stop disc radius
STOP_MARGIN = (3.8, 5.2) # free track after the first / before the last station
CORNER_MARGIN = 3.0      # no stop this close to a bend
STOP_SPACING_MAX = 9.0   # stops on a long stretch are grouped, not spread out
LABEL_GAP = 2.3          # stop centre → label
CHEVRON = 1.9            # direction chevron size
STRIPE = (1.3, 1.3)      # partial step: full-colour stripe / 40 %-tint gap along the line
PARTIAL_ALPHA = 0.4      # opacity of the tinted base of a partial step


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
            "off": "",          # no step is switched off by genome type in this cycle
            "reversed": bool(st.get("reversible")
                             and (ctx or {}).get("dsr_direction") == "oxidative"),
            "lines": [[_token(s, codes) for s in line] for line in st["lines"]],
        }
    return out


# ───────────────────────────── geometry ──────────────────────────────────────

def _pt(p):
    """A waypoint: a station id or an (x, y) pair."""
    return (NODES[p][1], NODES[p][2]) if isinstance(p, str) else tuple(p)


def _offset_polyline(pts, d):
    """Shift a polyline by d to the left of its direction of travel (y up);
    corners are joined by the intersection of the shifted segments."""
    if not d or len(pts) < 2:
        return list(pts)
    segs = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        L = math.hypot(x1 - x0, y1 - y0) or 1.0
        nx, ny = -(y1 - y0) / L, (x1 - x0) / L
        segs.append(((x0 + nx * d, y0 + ny * d), (x1 + nx * d, y1 + ny * d)))
    out = [segs[0][0]]
    for (a0, a1), (b0, b1) in zip(segs, segs[1:]):
        ax, ay = a1[0] - a0[0], a1[1] - a0[1]
        bx, by = b1[0] - b0[0], b1[1] - b0[1]
        den = ax * by - ay * bx
        if abs(den) < 1e-9:
            out.append(a1)
            continue
        t = ((b0[0] - a0[0]) * by - (b0[1] - a0[1]) * bx) / den
        out.append((a0[0] + ax * t, a0[1] + ay * t))
    out.append(segs[-1][1])
    return out


def _along(pts, dist):
    """Point, unit tangent and stretch index at a distance along the polyline."""
    left = dist
    for k, ((x0, y0), (x1, y1)) in enumerate(zip(pts, pts[1:])):
        L = math.hypot(x1 - x0, y1 - y0)
        if left <= L or (x1, y1) == pts[-1]:
            t = min(max(left / L, 0.0), 1.0) if L else 0.0
            return (x0 + (x1 - x0) * t, y0 + (y1 - y0) * t), \
                (((x1 - x0) / L, (y1 - y0) / L) if L else (1.0, 0.0)), k
        left -= L
    return pts[-1], (1.0, 0.0), len(pts) - 2


def _length(pts) -> float:
    return sum(math.hypot(x1 - x0, y1 - y0) for (x0, y0), (x1, y1) in zip(pts, pts[1:]))


def _stop_distances(pts, n: int, shift: float = 0.0,
                    spacing: float = STOP_SPACING_MAX) -> list[float]:
    """Distances along the polyline of n evenly spread stops, keeping clear of
    the two stations and of every bend; a long stretch groups its stops around
    its middle instead of spreading them out."""
    m0, m1 = STOP_MARGIN
    usable, at = [], 0.0
    nseg = len(pts) - 1
    for k, ((x0, y0), (x1, y1)) in enumerate(zip(pts, pts[1:])):
        L = math.hypot(x1 - x0, y1 - y0)
        a = at + (m0 if k == 0 else CORNER_MARGIN)
        b = at + L - (m1 if k == nseg - 1 else CORNER_MARGIN)
        if b > a:
            usable.append((a, b))
        at += L
    U = sum(b - a for a, b in usable) or 1.0
    span = min(U, n * spacing)
    lead = min(max((U - span) / 2 + shift, 0.0), U - span)
    out = []
    for i in range(n):
        u = lead + span * (i + 0.5) / n
        for a, b in usable:
            if u <= b - a + 1e-9:
                out.append(a + u)
                break
            u -= b - a
        else:
            out.append(usable[-1][1] if usable else L / 2)
    return out


def stripes(pts, width: float = LINE_W) -> list[list[tuple[float, float]]]:
    """45° hatching of a partial step: parallelograms of full colour across the
    line, drawn over the tinted base. Each straight stretch is hatched on its
    own, keeping one line width clear of the ends and bends."""
    s_len, gap = STRIPE
    quads = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        L = math.hypot(x1 - x0, y1 - y0)
        if L < 2 * width + s_len:
            continue
        tx, ty = (x1 - x0) / L, (y1 - y0) / L
        nx, ny = -ty, tx
        h = width / 2
        at = lambda a, b: (x0 + tx * a + nx * b, y0 + ty * a + ny * b)   # noqa: E731
        c = width + s_len / 2
        while c + s_len / 2 + h <= L - width:
            quads.append([at(c - s_len / 2 - h, -h), at(c + s_len / 2 - h, -h),
                          at(c + s_len / 2 + h, h), at(c - s_len / 2 + h, h)])
            c += s_len + gap
    return quads


def _orient(tangent) -> str:
    tx, ty = tangent
    if abs(ty) > 0.9:
        return "v"
    if abs(tx) > 0.9:
        return "h"
    return "d"


def step_geometry(st) -> dict:
    """Polyline, stops (one per gene token, in `lines` order) and the chevron
    positions of one step. Stop labels: on a vertical or diagonal stretch the
    text is horizontal beside the stop; on a horizontal stretch it is written
    at 45° so dense stops do not collide (the metro-map convention), unless
    the step is `flat` (a lone stop between parallel strands: horizontal text
    centred over or under the stop)."""
    pts = _offset_polyline([_pt(p) for p in st["path"]], st.get("offset", 0.0))
    L = _length(pts)
    tokens = [g for line in st["lines"] for g in line]
    n = len(tokens)
    stops = []
    sides = st.get("side", 1)
    for i, d in enumerate(_stop_distances(pts, n, st.get("stop_shift", 0.0),
                                          st.get("stop_spacing", STOP_SPACING_MAX))):
        (x, y), (tx, ty), k = _along(pts, d)
        side = sides[k] if isinstance(sides, (list, tuple)) else sides
        side = side or (1 if i % 2 == 0 else -1)
        nx, ny = -ty * side, tx * side          # left normal × side
        o = _orient((tx, ty))
        if o == "h" and st.get("flat"):
            up = 1 if ny > 0 else -1
            stops.append({"x": x, "y": y, "lx": x, "ly": y + LABEL_GAP * up, "rot": 0,
                          "ha": "center", "va": "bottom" if up > 0 else "top"})
        elif o == "h":
            # 45° label: starts above the stop reading up-right, or ends below
            # it reading up to the stop.
            up = 1 if ny > 0 else -1
            lx, ly = x + 0.55 * up, y + LABEL_GAP * up
            stops.append({"x": x, "y": y, "lx": lx, "ly": ly, "rot": 45,
                          "ha": "left" if up > 0 else "right", "va": "center"})
        else:
            lx, ly = x + nx * LABEL_GAP, y + ny * LABEL_GAP
            ha = "left" if nx > 0.3 else "right" if nx < -0.3 else "center"
            va = "center" if o == "v" else ("bottom" if ny > 0 else "top")
            stops.append({"x": x, "y": y, "lx": lx, "ly": ly, "rot": 0,
                          "ha": ha, "va": va})

    # Chevron tips: just before each end station, pointing along the travel.
    def chevron(at_end: bool):
        d = L - STATION_R - 1.0 if at_end else STATION_R + 1.0
        (x, y), (tx, ty), _ = _along(pts, d)
        if not at_end:
            tx, ty = -tx, -ty
        return {"x": x, "y": y, "tx": tx, "ty": ty}
    feeders = []
    for fp in st.get("feeders", []):
        fpts = [_pt(p) for p in fp]
        (fx, fy), (ftx, fty), _ = _along(fpts, _length(fpts) - STATION_R - 1.0)
        feeders.append({"path": fpts, "stripes": stripes(fpts),
                        "chevron": {"x": fx, "y": fy, "tx": ftx, "ty": fty}})
    return {"path": pts, "stops": stops, "chevron": chevron(True),
            "chevron_rev": chevron(False), "stripes": stripes(pts),
            "feeders": feeders}


def layout() -> dict:
    """Everything a renderer needs, JSON-serialisable."""
    r = lambda v: round(v, 2)                                 # noqa: E731
    steps = []
    for st in STEPS:
        g = step_geometry(st)
        steps.append({
            "id": st["id"], "name": st["name"], "pathway": st["pathway"],
            "path": [[r(x), r(y)] for x, y in g["path"]],
            "stops": [{k: (r(v) if isinstance(v, float) else v) for k, v in s.items()}
                      for s in g["stops"]],
            "chevron": {k: r(v) for k, v in g["chevron"].items()},
            "chevron_rev": {k: r(v) for k, v in g["chevron_rev"].items()},
            "stripes": [[[r(x), r(y)] for x, y in q] for q in g["stripes"]],
            "feeders": [{"path": [[r(x), r(y)] for x, y in fd["path"]],
                         "stripes": [[[r(x), r(y)] for x, y in q] for q in fd["stripes"]],
                         "chevron": {k: r(v) for k, v in fd["chevron"].items()}}
                        for fd in g["feeders"]],
        })
    return {
        "xlim": list(XLIM), "ylim": list(YLIM),
        "wide": (XLIM[1] - XLIM[0]) > 130,
        "geom": {"line_w": LINE_W, "station_r": STATION_R, "stop_r": STOP_R,
                 "chevron": CHEVRON, "partial_alpha": PARTIAL_ALPHA},
        "nodes": [{"id": k, "label": v[0], "x": v[1], "y": v[2],
                   "ldx": v[3][0], "ldy": v[3][1], "ha": v[3][2], "va": v[3][3]}
                  for k, v in NODES.items()],
        "steps": steps,
    }
