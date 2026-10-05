"""Everything specific to the SULFUR cycle that the figures and report need:
pathway order, palette, labels and output file names. The nitrogen sister
pipeline has its own _domain.py and _cycle_model.py; every other figure /
report script is identical between the two.

Palette: the pathway hues are the first six slots of the same validated
palette the nitrogen pipeline uses, assigned in CAT_ORDER (the order the
pathway blocks sit next to each other; adjacent CVD ΔE >= 8, normal-vision
ΔE >= 15, both modes). Three light-mode hues are below 3:1 on white, so colour
never carries identity alone: every coloured mark has a text label or a block
header.
"""

CYCLE_LETTER = "S"
CYCLE_NAME = "sulfur"
CALLS_TSV = "scycle_calls.tsv"
LOCI_TSV = "scycle_loci.tsv"
TOOL = "scycle"
REPORT_HTML = "scycle_report.html"

CAT_ORDER = [
    "dissimilatory_sulfate_reduction", "sulfur_oxidation",
    "assimilatory_sulfate_reduction", "thiosulfate_polysulfide",
    "organic_sulfur_dmsp", "sulfonate_taurine",
]

PATHWAY_COLOR = {
    "dissimilatory_sulfate_reduction": "#2a78d6",   # blue
    "sulfur_oxidation":                "#eb6834",   # orange
    "assimilatory_sulfate_reduction":  "#1baf7a",   # aqua
    "thiosulfate_polysulfide":         "#eda100",   # yellow
    "organic_sulfur_dmsp":             "#e87ba4",   # magenta
    "sulfonate_taurine":               "#008300",   # green
}
# Same hues stepped for a dark surface (used by report.html only).
PATHWAY_COLOR_DARK = {
    "dissimilatory_sulfate_reduction": "#3987e5",
    "sulfur_oxidation":                "#d95926",
    "assimilatory_sulfate_reduction":  "#199e70",
    "thiosulfate_polysulfide":         "#c98500",
    "organic_sulfur_dmsp":             "#d55181",
    "sulfonate_taurine":               "#008300",
}
PATHWAY_LABEL = {
    "dissimilatory_sulfate_reduction": "Dissimilatory sulfate reduction / reverse Dsr",
    "sulfur_oxidation":                "Sulfur oxidation",
    "assimilatory_sulfate_reduction":  "Assimilatory sulfate reduction",
    "thiosulfate_polysulfide":         "Thiosulfate / tetrathionate",
    "organic_sulfur_dmsp":             "Organic S (DMSP / DMS)",
    "sulfonate_taurine":               "Sulfonate / taurine",
}
# Shorter form for the angled block headers of the overview grid.
PATHWAY_SHORT = dict(
    PATHWAY_LABEL,
    dissimilatory_sulfate_reduction="Dissimilatory SO₄²⁻ reduction",
    assimilatory_sulfate_reduction="Assimilatory SO₄²⁻ reduction",
)

# Tidy-ups applied to complex / module ids when they are shown as labels.
PRETTY_REPLACE = [("aps ", "APS "), ("dsr", "Dsr"), ("sox", "Sox"),
                  ("qmo", "Qmo"), ("soe ", "Soe "), ("atp ", "ATP "),
                  ("dmsp", "DMSP"), (" dh", " dehydrogenase"),
                  ("assim ", "assimilatory ")]

# Legend text for a module that is ruled out by an exclusion rule.
RULED_OUT_LABEL = "ruled out — dsrAB runs in the other direction"
