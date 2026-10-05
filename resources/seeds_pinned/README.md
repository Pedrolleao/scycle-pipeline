# Pinned BLAST seeds

`unstable_refs.fasta` (38 sequences) and `blast_gated_refs.fasta` (56 sequences):
the UniProt reference sequences behind the validated DIAMOND / BLAST databases, as
fetched up to 2026-05-27. Headers are `>{target_id}||{accession}` followed by the
UniProt header.

**Why they are in the repository.** UniProt entries are revised and deleted, so a
database rebuilt from UniProt is not guaranteed to be the validated one.
`workflow/scripts/build_blast_db.py` copies these files into `resources/blast_db/` and
indexes them. With `--upstream`, or when `config/targets.yaml` lists a curated
accession that is not in these files, it fetches all seeds from UniProt instead —
after adding a seed, rebuild that way, check the result, and replace the two files here.

Sequences are from UniProtKB (CC BY 4.0).
