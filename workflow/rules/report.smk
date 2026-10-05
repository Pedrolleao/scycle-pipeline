# report.smk — bridges per-sample calls to:
#   1. per-sample complex/synergy completeness + gap_analysis.txt
#      + S-cycle map + (nucleotide input) gene-neighbourhood maps
#   2. cross-sample matrix.tsv + heatmap + focused figures + scycle_report.html


rule complex_completeness:
    input:
        calls=RESULTS / "{sample}" / "calls" / "scycle_calls.tsv",
        targets=config["paths"]["targets_yaml"],
    output:
        complexes=RESULTS / "{sample}" / "calls" / "complex_completeness.tsv",
        synergies=RESULTS / "{sample}" / "calls" / "synergy_completeness.tsv",
    shell:
        r"""
        python workflow/scripts/compute_complex_completeness.py \
            --calls {input.calls} --targets {input.targets} \
            --complexes-out {output.complexes} \
            --synergies-out {output.synergies}
        """


rule gap_analysis:
    input:
        calls=RESULTS / "{sample}" / "calls" / "scycle_calls.tsv",
        complexes=RESULTS / "{sample}" / "calls" / "complex_completeness.tsv",
        synergies=RESULTS / "{sample}" / "calls" / "synergy_completeness.tsv",
        targets=config["paths"]["targets_yaml"],
    output:
        txt=RESULTS / "{sample}" / "report" / "gap_analysis.txt",
    params:
        sample="{sample}",
        mode=lambda wc: "read" if sample_is_read_mode(wc.sample) else "protein",
    shell:
        r"""
        python workflow/scripts/make_gap_analysis.py \
            --sample {params.sample} --mode {params.mode} \
            --calls {input.calls} \
            --complexes {input.complexes} \
            --synergies {input.synergies} \
            --targets {input.targets} \
            --out {output.txt}
        """


rule cross_sample_report:
    input:
        calls=[str(RESULTS / s / "calls" / "scycle_calls.tsv") for s in SAMPLES],
        complexes=[str(RESULTS / s / "calls" / "complex_completeness.tsv") for s in SAMPLES],
        synergies=[str(RESULTS / s / "calls" / "synergy_completeness.tsv") for s in SAMPLES],
        targets=config["paths"]["targets_yaml"],
    output:
        matrix=RESULTS / "scycle_matrix.tsv",
        heatmap_svg=RESULTS / "scycle_heatmap.svg",
        heatmap_png=RESULTS / "scycle_heatmap.png",
    params:
        samples=",".join(SAMPLES.keys()),
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_cross_sample_report.py \
            --samples {params.samples} \
            --results-dir {params.results_dir} \
            --targets {input.targets} \
            --matrix-out {output.matrix} \
            --heatmap-out {output.heatmap_svg}
        python workflow/scripts/make_cross_sample_report.py \
            --samples {params.samples} \
            --results-dir {params.results_dir} \
            --targets {input.targets} \
            --matrix-out {output.matrix} \
            --heatmap-out {output.heatmap_png}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Focused per-pathway / complex / synergy / S-cycle-map figures.
# Replace the legacy dense `scycle_heatmap.svg` for daily reading — the
# legacy heatmap is still produced above as a quick-glance overview.
#
# PATHWAY_CATEGORIES is defined in the top-level Snakefile.
# ─────────────────────────────────────────────────────────────────────────────

FIGURES_DIR = RESULTS / "figures"


rule pathway_heatmap:
    """One categorical heatmap per pathway category. Co-renders SVG (vector,
    paper-grade) and PNG (300 DPI, quick preview)."""
    input:
        matrix=RESULTS / "scycle_matrix.tsv",
        targets=config["paths"]["targets_yaml"],
    output:
        svg=FIGURES_DIR / "pathway_{pathway}.svg",
        png=FIGURES_DIR / "pathway_{pathway}.png",
    wildcard_constraints:
        pathway="|".join(PATHWAY_CATEGORIES),
    shell:
        r"""
        python workflow/scripts/make_pathway_heatmap.py \
            --matrix {input.matrix} --targets {input.targets} \
            --pathway {wildcards.pathway} --out {output.svg}
        python workflow/scripts/make_pathway_heatmap.py \
            --matrix {input.matrix} --targets {input.targets} \
            --pathway {wildcards.pathway} --out {output.png}
        """


rule complex_heatmap_figure:
    """Per-sample obligatory-complex completeness panel. SVG + 300 DPI PNG."""
    input:
        matrix=RESULTS / "scycle_matrix.tsv",        # ensures samples present
        complexes=[str(RESULTS / s / "calls" / "complex_completeness.tsv") for s in SAMPLES],
        targets=config["paths"]["targets_yaml"],
    output:
        svg=FIGURES_DIR / "complexes.svg",
        png=FIGURES_DIR / "complexes.png",
    params:
        samples=",".join(SAMPLES.keys()),
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_complex_heatmap.py \
            --results-dir {params.results_dir} --targets {input.targets} \
            --samples {params.samples} --out {output.svg}
        python workflow/scripts/make_complex_heatmap.py \
            --results-dir {params.results_dir} --targets {input.targets} \
            --samples {params.samples} --out {output.png}
        """


rule synergy_heatmap_figure:
    """Per-sample predicted-synergy completeness panel. SVG + 300 DPI PNG."""
    input:
        matrix=RESULTS / "scycle_matrix.tsv",
        synergies=[str(RESULTS / s / "calls" / "synergy_completeness.tsv") for s in SAMPLES],
        targets=config["paths"]["targets_yaml"],
    output:
        svg=FIGURES_DIR / "synergies.svg",
        png=FIGURES_DIR / "synergies.png",
    params:
        samples=",".join(SAMPLES.keys()),
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_synergy_heatmap.py \
            --results-dir {params.results_dir} --targets {input.targets} \
            --samples {params.samples} --out {output.svg}
        python workflow/scripts/make_synergy_heatmap.py \
            --results-dir {params.results_dir} --targets {input.targets} \
            --samples {params.samples} --out {output.png}
        """


rule scycle_maps_figure:
    """Small multiples: every genome's calls drawn on the sulfur cycle
    (replaces the former `applications` panel). SVG + 300 DPI PNG, one render."""
    input:
        calls=[str(RESULTS / s / "calls" / "scycle_calls.tsv") for s in SAMPLES],
        synergies=[str(RESULTS / s / "calls" / "synergy_completeness.tsv") for s in SAMPLES],
    output:
        svg=FIGURES_DIR / "scycle_maps.svg",
        png=FIGURES_DIR / "scycle_maps.png",
    params:
        samples=",".join(SAMPLES.keys()),
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_cycle_map.py \
            --results-dir {params.results_dir} --samples {params.samples} \
            --out {output.svg} {output.png}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Per-sample figures: the genome's N-cycle map, and (nucleotide / MAG input
# only — it needs gene coordinates) the gene-neighbourhood maps.
# ─────────────────────────────────────────────────────────────────────────────

rule scycle_map:
    """One genome's calls drawn on the sulfur cycle."""
    input:
        calls=RESULTS / "{sample}" / "calls" / "scycle_calls.tsv",
        synergies=RESULTS / "{sample}" / "calls" / "synergy_completeness.tsv",
    output:
        svg=RESULTS / "{sample}" / "report" / "scycle_map.svg",
        png=RESULTS / "{sample}" / "report" / "scycle_map.png",
    params:
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_cycle_map.py \
            --results-dir {params.results_dir} --sample {wildcards.sample} \
            --out {output.svg} {output.png}
        """


rule locus_maps:
    """Cluster the called genes into loci and draw their neighbourhoods."""
    input:
        calls=RESULTS / "{sample}" / "calls" / "scycle_calls.tsv",
        faa=RESULTS / "{sample}" / "prodigal" / "{sample}.faa",
        gff=RESULTS / "{sample}" / "prodigal" / "{sample}.gff",
        targets=config["paths"]["targets_yaml"],
    output:
        tsv=RESULTS / "{sample}" / "calls" / "scycle_loci.tsv",
        svg=RESULTS / "{sample}" / "report" / "loci.svg",
        png=RESULTS / "{sample}" / "report" / "loci.png",
    params:
        results_dir=str(RESULTS),
    wildcard_constraints:
        sample=NUCLEOTIDE_SAMPLES_RE,
    shell:
        r"""
        python workflow/scripts/make_locus_maps.py \
            --sample {wildcards.sample} --results-dir {params.results_dir} \
            --targets {input.targets} --tsv {output.tsv} \
            --out {output.svg} {output.png}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Interactive report: one self-contained HTML file for the whole run.
# ─────────────────────────────────────────────────────────────────────────────

rule html_report:
    input:
        calls=[str(RESULTS / s / "calls" / "scycle_calls.tsv") for s in SAMPLES],
        complexes=[str(RESULTS / s / "calls" / "complex_completeness.tsv") for s in SAMPLES],
        synergies=[str(RESULTS / s / "calls" / "synergy_completeness.tsv") for s in SAMPLES],
        coords=[str(RESULTS / s / "prodigal" / f"{s}.{ext}")
                for s in NUCLEOTIDE_SAMPLES for ext in ("faa", "gff")],
        targets=config["paths"]["targets_yaml"],
        template="workflow/scripts/report_template.html",
    output:
        html=RESULTS / "scycle_report.html",
    params:
        samples=",".join(SAMPLES.keys()),
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_html_report.py \
            --samples {params.samples} --results-dir {params.results_dir} \
            --targets {input.targets} --out {output.html}
        """
