#!/usr/bin/env python3
"""Query UniProt for canonical DsrA / DsrB per target organism, with length + name
filters, so we get REAL full-length subunits (not fragments or fusion proteins).

For each (organism, gene, type, subtype) we run a UniProt search restricted by
protein/gene name and sequence-length window, prefer reviewed, take best hit.
Writes a candidate TSV we then materialize into faa/tsv.
"""
import json, sys, time, urllib.request, urllib.parse, os

OUTDIR = os.path.dirname(os.path.abspath(__file__))
# (organism query, gene dsrA|dsrB, type, subtype, label_org)
# dsrA target len ~ 380-470 ; dsrB ~ 340-440
TARGETS = [
    # --- reductive bacterial ---
    ("Desulfovibrio vulgaris", "reductive", "reductive_bacterial"),
    ("Desulfovibrio desulfuricans", "reductive", "reductive_bacterial"),
    ("Desulfovibrio alaskensis", "reductive", "reductive_bacterial"),
    ("Desulfovibrio gigas", "reductive", "reductive_bacterial"),
    ("Desulfomicrobium baculatum", "reductive", "reductive_bacterial"),
    ("Desulfobacter", "reductive", "reductive_bacterial"),
    ("Desulfobacterium autotrophicum", "reductive", "reductive_bacterial"),
    ("Desulfobulbus", "reductive", "reductive_bacterial"),
    ("Desulfococcus oleovorans", "reductive", "reductive_bacterial"),
    ("Desulfotomaculum", "reductive", "reductive_bacterial"),
    ("Desulfitobacterium hafniense", "reductive", "reductive_bacterial"),
    ("Desulfosporosinus", "reductive", "reductive_bacterial"),
    ("Thermodesulfovibrio yellowstonii", "reductive", "reductive_bacterial"),
    ("Thermodesulfobacterium", "reductive", "reductive_bacterial"),
    ("Desulfarculus baarsii", "reductive", "reductive_bacterial"),
    ("Desulfurivibrio alkaliphilus", "reductive", "reductive_bacterial"),
    ("Desulfobacca acetoxidans", "reductive", "reductive_bacterial"),
    ("Desulfotignum phosphitoxidans", "reductive", "reductive_bacterial"),
    ("Candidatus Desulforudis audaxviator", "reductive", "reductive_bacterial"),
    # --- reductive archaeal ---
    ("Archaeoglobus fulgidus", "reductive", "reductive_archaeal"),
    ("Archaeoglobus profundus", "reductive", "reductive_archaeal"),
    ("Archaeoglobus veneficus", "reductive", "reductive_archaeal"),
    ("Archaeoglobus sulfaticallidus", "reductive", "reductive_archaeal"),
    # --- oxidative bacterial (reverse Dsr) ---
    ("Allochromatium vinosum", "oxidative", "oxidative_bacterial"),
    ("Chlorobaculum tepidum", "oxidative", "oxidative_bacterial"),
    ("Chlorobium limicola", "oxidative", "oxidative_bacterial"),
    ("Chlorobium phaeobacteroides", "oxidative", "oxidative_bacterial"),
    ("Thiobacillus denitrificans", "oxidative", "oxidative_bacterial"),
    ("Thioalkalivibrio", "oxidative", "oxidative_bacterial"),
    ("Thioflavicoccus mobilis", "oxidative", "oxidative_bacterial"),
    ("Thiocapsa", "oxidative", "oxidative_bacterial"),
    ("Sideroxydans lithotrophicus", "oxidative", "oxidative_bacterial"),
    ("Magnetococcus marinus", "oxidative", "oxidative_bacterial"),
    ("Beggiatoa", "oxidative", "oxidative_bacterial"),
    ("Thiothrix", "oxidative", "oxidative_bacterial"),
    ("Rhodovulum sulfidophilum", "oxidative", "oxidative_bacterial"),
    ("Halorhodospira", "oxidative", "oxidative_bacterial"),
]

LEN = {"dsrA": (360, 520), "dsrB": (320, 470)}

def uquery(q, lo, hi):
    fields = "accession,id,gene_names,protein_name,organism_name,length,reviewed,sequence"
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urllib.parse.urlencode({
               "query": q,
               "fields": fields,
               "format": "json",
               "size": "25",
           }))
    req = urllib.request.Request(url, headers={"User-Agent": "scycle-curation/1.0"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            sys.stderr.write(f"retry {attempt} {q[:40]}: {e}\n")
            time.sleep(2)
    return {"results": []}

def best(results, gene, lo, hi):
    """pick best: reviewed first, length in window, name matches gene."""
    cands = []
    for r in results.get("results", []):
        L = r.get("sequence", {}).get("length", 0)
        seq = r.get("sequence", {}).get("value", "")
        pname = (r.get("proteinDescription", {})
                  .get("recommendedName", {}).get("fullName", {}).get("value", "")) or ""
        # gather submitted names too
        allnames = pname.lower()
        gn = ""
        for g in r.get("genes", []):
            gn += (g.get("geneName", {}).get("value", "") + " ").lower()
        rev = r.get("entryType", "").lower().startswith("uniprotkb reviewed") or "reviewed" in r.get("entryType","").lower()
        acc = r.get("primaryAccession", "")
        org = r.get("organism", {}).get("scientificName", "")
        # subunit match
        sub = "alpha" if gene == "dsrA" else "beta"
        subok = (sub in allnames) or (gene.lower() in gn) or (gene.lower() in allnames)
        sulfite = "sulfite reductase" in allnames or "dsr" in gn or "dsr" in allnames
        if not (lo <= L <= hi):
            continue
        if not seq:
            continue
        score = 0
        if rev: score += 10
        if sulfite: score += 5
        if subok: score += 3
        cands.append((score, acc, org, L, seq, pname, gn.strip(), rev))
    cands.sort(reverse=True)
    return cands

def main():
    rows = []  # acc, gene, type, subtype, org, source, notes, seq
    seen_acc = set()
    log = []
    for orgq, typ, subtype in TARGETS:
        for gene in ("dsrA", "dsrB"):
            lo, hi = LEN[gene]
            sub = "alpha" if gene == "dsrA" else "beta"
            q = (f'(organism_name:"{orgq}") AND '
                 f'(protein_name:"dissimilatory sulfite reductase" OR gene:{gene} OR gene:dsrAB) '
                 f'AND (protein_name:{sub} OR gene:{gene})')
            res = uquery(q, lo, hi)
            cands = best(res, gene, lo, hi)
            time.sleep(0.4)
            if not cands:
                # fallback broader
                q2 = f'(organism_name:"{orgq}") AND (gene:{gene})'
                res = uquery(q2, lo, hi)
                cands = best(res, gene, lo, hi)
                time.sleep(0.4)
            picked = None
            for c in cands:
                if c[1] not in seen_acc:
                    picked = c
                    break
            if picked is None:
                log.append(f"NONE {orgq} {gene}")
                continue
            _, acc, org, L, seq, pname, gn, rev = picked
            seen_acc.add(acc)
            note = f"len={L};{'reviewed' if rev else 'unreviewed'};{pname[:40] or gn}"
            rows.append((acc, gene, typ, subtype, org, "UniProt", note, seq))
            log.append(f"OK {orgq} {gene} -> {acc} L={L}")
    # write candidate intermediate
    with open(os.path.join(OUTDIR, "_dsrab_candidates.tsv"), "w") as f:
        f.write("accession\tgene\ttype\tsubtype\torganism\tsource\tnotes\tlen\n")
        for acc, gene, typ, subtype, org, src, note, seq in rows:
            f.write(f"{acc}\t{gene}\t{typ}\t{subtype}\t{org}\t{src}\t{note}\t{len(seq)}\n")
    # store seqs
    with open(os.path.join(OUTDIR, "_dsrab_seqs.json"), "w") as f:
        json.dump([(acc, gene, typ, subtype, org, src, note, seq)
                   for acc, gene, typ, subtype, org, src, note, seq in rows], f)
    with open(os.path.join(OUTDIR, "_dsrab_query_log.txt"), "w") as f:
        f.write("\n".join(log))
    print(f"rows={len(rows)}")
    from collections import Counter
    print(Counter((g, t) for _, g, t, _, _, _, _, _ in rows))

if __name__ == "__main__":
    main()
