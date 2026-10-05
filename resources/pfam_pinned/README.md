# Pinned Pfam profiles

PF04358 and PF08679, the Pfam profile HMM(s) used as fallback signature for the targets of
`config/targets.yaml` that have no KO, as fetched from the InterPro API
(`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<id>/?annotation=hmm`) on 2026-05-26
and built into the validated database. `workflow/scripts/build_hmm_db.py` uses these
files; with `--upstream`, or for a Pfam id that is not here, it fetches the current
profile from InterPro. Pfam is released under CC0.
