# Pinned KOfam snapshot

The 58 KOfam profile HMMs (`profiles/K*.hmm`) and the matching rows of `ko_list`
(`ko_list.tsv`) for the KOs named in `config/targets.yaml`, taken unchanged from the
KOfam release the tool was validated on.

| | |
|---|---|
| source | `https://www.genome.jp/ftp/db/kofam/profiles.tar.gz` and `ko_list.gz` |
| fetched | 2026-05-24 |
| SHA-256 of that `profiles.tar.gz` | `b03d20b96254d0f04652102ea538087ac36adadf0f807ee7278f8637810ca15a` |

**Why it is in the repository.** genome.jp serves KOfam from one rolling URL and keeps
no old releases. The release of 2026-09-29 gives a different threshold for 49 of these
58 KOs (K16937: 70.07 → 215.57), so a database built from a fresh download is not the validated one.
`workflow/scripts/build_hmm_db.py` builds from this directory by default;
`--upstream` downloads the current release instead.

**Changing the KO set.** A KO added to `config/targets.yaml` is not here, and the build
then takes *all* KOs from the download (two releases are never mixed). To extend the
snapshot, add the profile and the `ko_list` row from the same pinned release.

KOfam is the work of the Kanehisa Laboratories (KEGG): Aramaki T. et al. (2020)
KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold.
*Bioinformatics* 36:2251–2252. Check the KEGG terms of use before redistributing these
files outside this repository.
