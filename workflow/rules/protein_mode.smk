# protein_mode.smk — Mode A: protein input (.faa) or nucleotide input
# (.fna) → prodigal → .faa → hmmscan + DIAMOND-blastp fallback → apply_rules.

SCRIPTS = "workflow/scripts"


def proteome_for(sample):
    """Return the .faa path. For protein samples it is the input itself;
    for nucleotide samples it is the Prodigal output."""
    if sample_kind(sample) == "protein":
        return sample_path(sample)
    return str(RESULTS / sample / "prodigal" / f"{sample}.faa")


rule prodigal:
    """Translate nucleotide input to predicted proteins."""
    input:
        lambda wc: sample_path(wc.sample),
    output:
        faa=RESULTS / "{sample}" / "prodigal" / "{sample}.faa",
        gff=RESULTS / "{sample}" / "prodigal" / "{sample}.gff",
    params:
        mode=lambda wc: SAMPLES[wc.sample].get("prodigal_mode", "single"),
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    log:
        RESULTS / "{sample}" / "prodigal" / "{sample}.log",
    shell:
        r"""
        prodigal -i {input} -a {output.faa} -f gff -o {output.gff} \
                 -p {params.mode} -q > {log} 2>&1
        """


rule hmmscan:
    """Scan a proteome against the concatenated ewaste HMM database.
    Output: domain-table (`--domtblout`) for downstream parsing."""
    input:
        faa=lambda wc: proteome_for(wc.sample),
        hmm=config["paths"]["hmm_db"],
    output:
        tbl=RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.tsv",
    params:
        evalue=config["thresholds"]["hmmer"]["evalue"],
        dom_e=config["thresholds"]["hmmer"]["dom_evalue"],
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    threads: 4
    log:
        RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.log",
    shell:
        r"""
        hmmscan --cpu {threads} -E {params.evalue} --domE {params.dom_e} \
                --domtblout {output.tbl} {input.hmm} {input.faa} > {log} 2>&1
        """


rule diamond_blastp_unstable:
    """DIAMOND-blastp the proteome against the unstable_refs database
    (9 weak/⚠️ Pfam targets)."""
    input:
        faa=lambda wc: proteome_for(wc.sample),
        db=Path(config["paths"]["blast_unstable"]).with_suffix(".dmnd"),
    output:
        tsv=RESULTS / "{sample}" / "blast" / "{sample}.unstable.tsv",
    params:
        evalue=config["thresholds"]["blast"]["evalue"],
        qcov=config["thresholds"]["blast"]["qcov_hsp_perc"],
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    threads: 4
    log:
        RESULTS / "{sample}" / "blast" / "{sample}.unstable.log",
    shell:
        r"""
        diamond blastp --quiet --threads {threads} --evalue {params.evalue} \
                       --query-cover {params.qcov} \
                       --db {input.db} --query {input.faa} \
                       --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                                   qstart qend sstart send evalue bitscore \
                       --out {output.tsv} 2> {log}
        """


rule diamond_blastp_gated:
    """DIAMOND-blastp the proteome against the blast_gated_refs database
    (no-Pfam targets + targets with `requires_blast_for_confirmation: true`).
    Permissive pident threshold applied downstream in apply_rules.py."""
    input:
        faa=lambda wc: proteome_for(wc.sample),
        db=Path(config["paths"]["blast_gated"]).with_suffix(".dmnd"),
    output:
        tsv=RESULTS / "{sample}" / "blast" / "{sample}.blast_gated.tsv",
    params:
        evalue=config["thresholds"]["blast"]["evalue"],
        qcov=config["thresholds"]["blast"]["qcov_hsp_perc"],
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    threads: 4
    log:
        RESULTS / "{sample}" / "blast" / "{sample}.blast_gated.log",
    shell:
        r"""
        diamond blastp --quiet --threads {threads} --evalue {params.evalue} \
                       --query-cover {params.qcov} \
                       --db {input.db} --query {input.faa} \
                       --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                                   qstart qend sstart send evalue bitscore \
                       --out {output.tsv} 2> {log}
        """


rule apply_rules_protein:
    """Combine hmmscan + DIAMOND-blastp evidence and assign per-target status."""
    input:
        hmm=RESULTS / "{sample}" / "hmm" / "{sample}.hmmscan.tsv",
        blast_u=RESULTS / "{sample}" / "blast" / "{sample}.unstable.tsv",
        blast_g=RESULTS / "{sample}" / "blast" / "{sample}.blast_gated.tsv",
        targets=config["paths"]["targets_yaml"],
        # The proteome .faa doubles as the gene-coordinate source: for nucleotide
        # samples it is the Prodigal output (coords in the headers → nxr/nar synteny
        # resolution); for protein samples it is the pre-called input (no coords).
        faa=lambda wc: proteome_for(wc.sample),
    output:
        calls=RESULTS / "{sample}" / "calls" / "scycle_calls.tsv",
    params:
        sample="{sample}",
        mode="protein",
    wildcard_constraints:
        sample=PROTEIN_SAMPLES_RE,
    log:
        RESULTS / "{sample}" / "calls" / "apply_rules.log",
    shell:
        r"""
        python {SCRIPTS}/apply_rules.py \
            --sample {params.sample} --mode {params.mode} \
            --hmm {input.hmm} \
            --blast-unstable {input.blast_u} \
            --blast-gated {input.blast_g} \
            --targets {input.targets} \
            --gene-coords {input.faa} \
            --out {output.calls} > {log} 2>&1
        """


