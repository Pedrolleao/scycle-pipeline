# report.smk — bridges per-sample calls to:
#   1. per-sample complex/synergy completeness + gap_analysis.txt
#   2. cross-sample matrix.tsv + heatmap.svg


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
        matrix=RESULTS / "multisample_matrix.tsv",
        heatmap_svg=RESULTS / "multisample_heatmap.svg",
        heatmap_png=RESULTS / "multisample_heatmap.png",
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
# Focused per-pathway / complex / synergy / applications figures.
# Replace the legacy dense `multisample_heatmap.svg` for daily reading — the
# legacy heatmap is still produced above as a quick-glance overview.
#
# PATHWAY_CATEGORIES is defined in the top-level Snakefile.
# ─────────────────────────────────────────────────────────────────────────────

FIGURES_DIR = RESULTS / "figures"


rule pathway_heatmap:
    """One categorical heatmap per pathway category. Co-renders SVG (vector,
    paper-grade) and PNG (300 DPI, quick preview)."""
    input:
        matrix=RESULTS / "multisample_matrix.tsv",
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
        matrix=RESULTS / "multisample_matrix.tsv",        # ensures samples present
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
        matrix=RESULTS / "multisample_matrix.tsv",
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


rule applications_panel:
    """Per-genome small-multiples breakdown of the dominant industrial
    application call (the winning complex + its subunit status strip).
    SVG + 300 DPI PNG."""
    input:
        matrix=RESULTS / "multisample_matrix.tsv",
        complexes=[str(RESULTS / s / "calls" / "complex_completeness.tsv") for s in SAMPLES],
        targets=config["paths"]["targets_yaml"],
    output:
        svg=FIGURES_DIR / "applications.svg",
        png=FIGURES_DIR / "applications.png",
    params:
        results_dir=str(RESULTS),
    shell:
        r"""
        python workflow/scripts/make_applications_panel.py \
            --matrix {input.matrix} --results-dir {params.results_dir} \
            --targets {input.targets} --out {output.svg}
        python workflow/scripts/make_applications_panel.py \
            --matrix {input.matrix} --results-dir {params.results_dir} \
            --targets {input.targets} --out {output.png}
        """
