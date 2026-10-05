# scycle-pipeline make targets.
#
# `make regression` is the accuracy gate (scores the 43-genome sp1_panel vs
# ground truth and FAILs below the floors in validation/test_regression.py).
# It is the real test for this pipeline.
#
# `make test_protein` is a lighter end-to-end SMOKE test: it runs the pipeline on
# the reference genomes downloaded by `make test_data` and confirms a non-empty
# results matrix — i.e. the pipeline runs start-to-finish. It does NOT check call
# accuracy; that is `make regression`.

TEST_DATA := test_data
REF_GENOMES := \
    $(TEST_DATA)/Aferrooxidans_ATCC23270.fna \
    $(TEST_DATA)/Mextorquens_AM1.fna \
    $(TEST_DATA)/Ecoli_K12.fna \
    $(TEST_DATA)/Soneidensis_MR1.fna \
    $(TEST_DATA)/Cmetallidurans_CH34.fna

NCBI_BASE := https://ftp.ncbi.nlm.nih.gov/genomes/all

URL_AFERR  := $(NCBI_BASE)/GCF/000/021/485/GCF_000021485.1_ASM2148v1/GCF_000021485.1_ASM2148v1_genomic.fna.gz
URL_MEXT   := $(NCBI_BASE)/GCF/000/022/685/GCF_000022685.1_ASM2268v1/GCF_000022685.1_ASM2268v1_genomic.fna.gz
URL_ECOLI  := $(NCBI_BASE)/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
URL_SHEW   := $(NCBI_BASE)/GCF/000/146/165/GCF_000146165.2_ASM14616v2/GCF_000146165.2_ASM14616v2_genomic.fna.gz
URL_CUPR   := $(NCBI_BASE)/GCF/000/196/015/GCF_000196015.1_ASM19601v1/GCF_000196015.1_ASM19601v1_genomic.fna.gz

PANEL := sp1_panel   # 43-genome validation panel (proteomes + 1 .fna, in the repository)

.PHONY: env dbs test_data test_protein test ground-truth regression regression-score benchmark scycdb clean clean_all

env:
	@mamba env create -f envs/scycle.yaml 2>/dev/null || conda env create -f envs/scycle.yaml

test_data: $(REF_GENOMES)

$(TEST_DATA)/Aferrooxidans_ATCC23270.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_AFERR) | gunzip > $@

$(TEST_DATA)/Mextorquens_AM1.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_MEXT) | gunzip > $@

$(TEST_DATA)/Ecoli_K12.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_ECOLI) | gunzip > $@

$(TEST_DATA)/Soneidensis_MR1.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_SHEW) | gunzip > $@

$(TEST_DATA)/Cmetallidurans_CH34.fna:
	mkdir -p $(TEST_DATA)
	curl -L $(URL_CUPR) | gunzip > $@

# End-to-end smoke test: run the pipeline on the reference genomes and confirm a
# non-empty results matrix (i.e. the pipeline runs start-to-finish). Accuracy is
# checked separately by `make regression`.
test_protein: test_data
	python run.py --input $(TEST_DATA) --mode protein --prodigal-mode single
	@python -c "import csv, sys, pathlib; \
m = pathlib.Path('results/multisample_matrix.tsv'); \
sys.exit('FAIL: results/multisample_matrix.tsv missing — pipeline did not finish') if not m.exists() else None; \
rows = list(csv.reader(m.open(), delimiter='\t')); \
sys.exit('FAIL: results matrix is empty') if len(rows) < 2 else None; \
print(f'OK: pipeline ran end-to-end — {len(rows)-1} data rows x {len(rows[0])-1} columns in results/multisample_matrix.tsv'); \
print('    (call accuracy is checked by: make regression)')"

# (Re)build the HMM and BLAST databases from config/targets.yaml and the pinned
# snapshots under resources/ (KOfam, Pfam, BLAST seeds) — no download.
dbs:
	python workflow/scripts/build_hmm_db.py --force
	python workflow/scripts/build_blast_db.py --force

# ── Accuracy regression gate (the real test for this pipeline) ────────────────
# Full gate: run the pipeline on the reference panel, score against the committed
# ground truths, and FAIL (non-zero exit) if accuracy dropped below the floors in
# validation/test_regression.py. Needs the databases (`make dbs`, or any run).
regression:
	python run.py --input $(PANEL) --skip-db-setup --prodigal-mode single --cores 8
	python validation/score_scycle.py
	python validation/test_regression.py

# Fast gate: re-score the EXISTING results/ matrix and check the floors (no
# pipeline run). Use after a scoring/ground-truth change when calls are current.
# Scores against the curated-FUNCTION GT (the default); GT_FILE=ground_truth.tsv for KEGG-KO
# (written to scycle_{metrics,confusion}.ground_truth.tsv, next to the default tables).
regression-score:
	python validation/score_scycle.py
	python validation/test_regression.py

test: regression

# Rebuild both ground truths: KEGG-derived (queries the KEGG REST API for the 43
# organisms unless validation/.kegg_cache/ holds them — KEGG's annotation moves, so
# the result can differ from the committed table), then curated-function.
ground-truth:
	python validation/build_ground_truth.py
	python validation/build_curated_function_gt.py

# ── Head-to-head benchmark vs raw KofamScan (curated-FUNCTION GT) ─────────────
# Builds the curated-function GT, then runs the genome-cluster-bootstrap benchmark
# (scycle vs raw KofamScan) on it. Needs the panel scored into results/ first
# (run `make regression` or the pipeline once). The KEGG-KO GT contrast is the
# GT_FILE=ground_truth.tsv variant documented in validation/benchmark/README.md.
benchmark:
	python validation/build_curated_function_gt.py
	cd validation/benchmark && GT_FILE=curated_function_gt.tsv python benchmark_stats.py --bootstrap 10000 \
	    --metabolic ../../../comparators/metabolic.tsv \
	    --dram ../../../comparators/dram.tsv \
	    --scycdb scycdb.tsv

# Regenerate the SCycDB domain-database comparator TSV (the sulfur analogue of
# ncycle's NCycDB). Runs SCycDB (Yu et al. 2020) on the 43-genome panel via DIAMOND;
# writes validation/benchmark/scycdb.tsv (consumed by `make benchmark` --scycdb).
# Needs the SCycDB DB staged under comparators/SCyc/data/ (364 MB FASTA from
# github.com/qichao1984/SCycDB; gitignored). See comparators/build_scycdb_tsv.py.
scycdb:
	python comparators/build_scycdb_tsv.py

clean:
	rm -rf results/*

clean_all: clean
	rm -rf resources/hmm/* resources/blast_db/* resources/nt_refs/*
