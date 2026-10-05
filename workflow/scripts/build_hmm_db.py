#!/usr/bin/env python3
"""
build_hmm_db.py — build the concatenated HMM database for the scycle-pipeline.

Primary backbone: **KOfam**. For every KEGG Orthology (KO) listed in
`targets[].ko` of config/targets.yaml, take that KO's profile HMM and its adaptive
score threshold from the pinned KOfam snapshot shipped in resources/kofam_pinned/
(the release the tool was validated on). If a KO is not in the snapshot, or with
`--upstream`, every KO is taken from the KOfam `profiles.tar.gz` / `ko_list`
downloaded from genome.jp instead — a rolling release, so thresholds may differ.
Pfam fallback profiles likewise come from resources/pfam_pinned/ when present.
This is the KO-primary detection tier the sulfur-cycle tool is built around.

Then:
  - splice in any `custom_hmm: true` target's targets/{tid}/{tid}.hmm
    (the homology-trap clade HMMs — added in the hardening phase; missing files
    are skipped with a warning, and apply_rules.py falls through to the KO tier);
  - for any target that has a `pfam:` but NO `ko:` (e.g. archaeal amoA, PF12942),
    fetch that Pfam HMM from InterPro as a fallback signature.

Outputs:
  resources/hmm/scycle_targets.hmm   (+ hmmpress index)
  resources/hmm/tc_cutoffs.tsv       (profile_id -> threshold; KO numbers,
                                       custom target ids, and Pfam ids)

Usage:  python workflow/scripts/build_hmm_db.py [--force] [--upstream]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

import yaml

from _common import ROOT, load_targets

HMM_DIR = ROOT / "resources" / "hmm"
CACHE = ROOT / "resources" / ".cache"
PINNED = ROOT / "resources" / "kofam_pinned"   # profiles/{KO}.hmm + ko_list.tsv
PFAM_PINNED = ROOT / "resources" / "pfam_pinned"   # {PFxxxxx}.hmm
CONCAT = HMM_DIR / "scycle_targets.hmm"
TC_TSV = HMM_DIR / "tc_cutoffs.tsv"
TARGETS_DIR = ROOT / "targets"

KOFAM_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
KOFAM_KO_LIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"
INTERPRO_HMM_URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pf}/?annotation=hmm"

# ── Pinned KOfam release (reproducibility) ────────────────────────────────────
# KEGG serves profiles.tar.gz from a rolling-latest URL with no embedded version,
# so we pin the exact release the validation campaign (curated_v3, F1 0.97) was
# built against by content hash. The build verifies the cached/downloaded tarball
# against this digest and warns loudly on a mismatch (= upstream rotated KOfam;
# scores may not reproduce). Bump both fields together after re-validating.
KOFAM_RELEASE = "2026-05-24 (genome.jp ftp latest at fetch time)"
KOFAM_PROFILES_SHA256 = "b03d20b96254d0f04652102ea538087ac36adadf0f807ee7278f8637810ca15a"


# ───────────────────────────── target inventory ──────────────────────────────

def unique_kos(targets: list[dict]) -> list[str]:
    seen, out = set(), []
    for t in targets:
        for ko in t.get("ko") or []:
            if ko not in seen:
                seen.add(ko)
                out.append(ko)
    return out


def pfam_only_ids(targets: list[dict]) -> list[str]:
    """Pfam IDs of targets that have NO KO — these need the Pfam fallback."""
    seen, out = set(), []
    for t in targets:
        if t.get("ko"):
            continue
        for pf in t.get("pfam") or []:
            if pf not in seen:
                seen.add(pf)
                out.append(pf)
    return out


# ───────────────────────────── KOfam snapshot ────────────────────────────────

def pinned_covers(kos: list[str]) -> bool:
    """True if resources/kofam_pinned/ holds a profile for every KO of targets.yaml.
    All or nothing: profiles of two KOfam releases are never mixed in one database."""
    return (PINNED / "ko_list.tsv").exists() and all(
        (PINNED / "profiles" / f"{ko}.hmm").exists() for ko in kos)


# ───────────────────────────── KOfam cache ───────────────────────────────────

def _sha256(path: Path, _buf: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_buf), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_kofam_release(profiles: Path) -> None:
    """Verify the profiles tarball against the pinned release digest and write a
    provenance stamp. A mismatch is a loud warning (not a hard fail): the user may
    be intentionally tracking a newer KOfam, but their scores then differ from the
    validated curated_v3 baseline."""
    digest = _sha256(profiles)
    stamp = CACHE / "kofam_release.txt"
    if digest == KOFAM_PROFILES_SHA256:
        print(f"[build_hmm_db] KOfam release verified — pinned {KOFAM_RELEASE}")
    else:
        print("[build_hmm_db] ! WARNING: KOfam profiles.tar.gz does NOT match the "
              "pinned release.", file=sys.stderr)
        print(f"               pinned sha256 {KOFAM_PROFILES_SHA256} ({KOFAM_RELEASE})",
              file=sys.stderr)
        print(f"               cached sha256 {digest}", file=sys.stderr)
        print("               Results may not reproduce the curated_v3 / F1 0.97 "
              "baseline. Re-validate and bump KOFAM_PROFILES_SHA256 if intended.",
              file=sys.stderr)
    stamp.write_text(
        f"pinned_release\t{KOFAM_RELEASE}\n"
        f"pinned_sha256\t{KOFAM_PROFILES_SHA256}\n"
        f"cached_sha256\t{digest}\n"
        f"profiles_url\t{KOFAM_PROFILES_URL}\n"
        f"matches_pin\t{digest == KOFAM_PROFILES_SHA256}\n")


def ensure_kofam_cache() -> tuple[Path, Path]:
    """Download ko_list + profiles.tar.gz into the cache once, verify the profiles
    tarball against the pinned KOfam release, and return their paths."""
    CACHE.mkdir(parents=True, exist_ok=True)
    ko_list = CACHE / "ko_list"
    ko_list_gz = CACHE / "ko_list.gz"
    profiles = CACHE / "profiles.tar.gz"
    if not ko_list.exists():
        if not ko_list_gz.exists():
            print(f"  downloading {KOFAM_KO_LIST_URL}…", file=sys.stderr)
            urllib.request.urlretrieve(KOFAM_KO_LIST_URL, ko_list_gz)
        with gzip.open(ko_list_gz, "rb") as src, open(ko_list, "wb") as dst:
            dst.write(src.read())
    if not profiles.exists():
        print(f"  downloading {KOFAM_PROFILES_URL} (~1.5 GB, one-time)…", file=sys.stderr)
        urllib.request.urlretrieve(KOFAM_PROFILES_URL, profiles)
    verify_kofam_release(profiles)
    return ko_list, profiles


def load_ko_thresholds(ko_list: Path) -> dict[str, tuple[float, str]]:
    """Return {KO: (threshold, score_type)} for KOs that carry a numeric threshold.
    KOs whose threshold is '-' (no profile threshold) are omitted — they then
    fall back to the hmmscan -E filter."""
    out: dict[str, tuple[float, str]] = {}
    with open(ko_list) as fh:
        next(fh, None)  # header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            ko, thr, stype = p[0], p[1], p[2]
            if thr in ("", "-"):
                continue
            try:
                out[ko] = (float(thr), stype)
            except ValueError:
                pass
    return out


def extract_ko_hmms(needed: set[str], profiles_tar: Path, outdir: Path) -> set[str]:
    """Extract profiles/{KO}.hmm for each needed KO from the KOfam tarball."""
    outdir.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    wanted = {f"profiles/{ko}.hmm" for ko in needed} | {f"{ko}.hmm" for ko in needed}
    with tarfile.open(profiles_tar, "r:gz") as tf:
        for m in tf:
            name = m.name.lstrip("./")
            if name in wanted:
                ko = Path(name).stem
                fobj = tf.extractfile(m)
                if fobj is None:
                    continue
                (outdir / f"{ko}.hmm").write_bytes(fobj.read())
                found.add(ko)
                if found >= needed:
                    break
    return found


# ───────────────────────────── Pfam fallback ─────────────────────────────────

def fetch_hmm_interpro(pf: str, out: Path) -> bool:
    url = INTERPRO_HMM_URL.format(pf=pf)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = r.read()
        try:
            text = gzip.decompress(data).decode()
        except OSError:
            text = data.decode(errors="replace")
        if "HMMER3" not in text.split("\n", 1)[0]:
            return False
        out.write_text(text)
        return True
    except Exception as e:
        print(f"  ! InterPro fetch failed for {pf}: {e}", file=sys.stderr)
        return False


def extract_tc(hmm_path: Path) -> str | None:
    """Return the TC (trusted cutoff) bitscore from an HMM header, if present."""
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith("TC "):
                return line.split()[1]
            if line.startswith("HMM "):
                break
    return None


# ───────────────────────────── main ──────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="Rebuild even if the concatenated DB exists")
    ap.add_argument("--upstream", action="store_true",
                    help="Ignore resources/kofam_pinned/ and pfam_pinned/ and download "
                         "the current KOfam release from genome.jp (~1.5 GB) and the "
                         "current Pfam profiles from InterPro")
    args = ap.parse_args()

    targets = load_targets()["targets"]
    kos = unique_kos(targets)
    pfams = pfam_only_ids(targets)
    # Optional per-target KO threshold override (targets.yaml `ko_tc:`) — used when
    # the KOfam ko_list threshold is mis-calibrated for an environmentally relevant
    # clade (e.g. nrfH K15876 too high for gammaproteobacterial NrfH).
    ko_override = {ko: float(t["ko_tc"]) for t in targets if t.get("ko_tc")
                   for ko in (t.get("ko") or [])}
    print(f"[build_hmm_db] {len(kos)} KOfam KOs + {len(pfams)} Pfam-only fallback profiles")

    HMM_DIR.mkdir(parents=True, exist_ok=True)
    if CONCAT.exists() and not args.force:
        print(f"[build_hmm_db] {CONCAT} exists; pass --force to rebuild")
        return

    if pinned_covers(kos) and not args.upstream:
        print(f"[build_hmm_db] using the pinned KOfam snapshot in {PINNED.relative_to(ROOT)}/ "
              f"({KOFAM_RELEASE})")
        ko_thr = load_ko_thresholds(PINNED / "ko_list.tsv")
        ko_hmm_dir = PINNED / "profiles"
        found = set(kos)
    else:
        if not args.upstream:
            print("[build_hmm_db] ! resources/kofam_pinned/ does not cover every KO of "
                  "targets.yaml — taking all of them from the KOfam download instead",
                  file=sys.stderr)
        ko_list, profiles_tar = ensure_kofam_cache()
        ko_thr = load_ko_thresholds(ko_list)
        ko_hmm_dir = CACHE / "ko_hmms"
        print(f"[build_hmm_db] extracting {len(kos)} KO profiles from {profiles_tar.name}…",
              flush=True)
        found = extract_ko_hmms(set(kos), profiles_tar, ko_hmm_dir)
    missing_ko = sorted(set(kos) - found)
    if missing_ko:
        print(f"[build_hmm_db] ! KO profiles not in tarball: {missing_ko}", file=sys.stderr)

    tc_rows = ["profile_id\tga\ttc"]
    with open(CONCAT, "w") as out:
        # 1. KOfam profiles.
        for ko in kos:
            hmm = ko_hmm_dir / f"{ko}.hmm"
            if not hmm.exists():
                continue
            out.write(hmm.read_text().rstrip("\n") + "\n")
            if ko in ko_override:
                tc_val = ko_override[ko]
            else:
                thr = ko_thr.get(ko)
                tc_val = thr[0] if thr else ""
            tc_rows.append(f"{ko}\t\t{tc_val}")

        # 2. Custom HMMs (homology-trap clade models; built in the hardening phase).
        custom_ids = [t["id"] for t in targets if t.get("custom_hmm")]
        for tid in custom_ids:
            chmm = TARGETS_DIR / tid / f"{tid}.hmm"
            if not chmm.exists():
                print(f"  · custom_hmm: true for '{tid}' but {chmm} not built yet "
                      "— KO tier will be used instead", file=sys.stderr)
                continue
            out.write(chmm.read_text().rstrip("\n") + "\n")
            tc_str = ""
            mf = chmm.parent / "manifest.yaml"
            if mf.exists():
                try:
                    m = yaml.safe_load(open(mf)) or {}
                    if m.get("tc_bitscore") is not None:
                        tc_str = f"{float(m['tc_bitscore']):.1f}"
                except Exception:
                    pass
            tc_rows.append(f"{tid}\t\t{tc_str}")
            print(f"  + spliced custom HMM: {tid} (TC={tc_str or 'unset'})")

        # 3. Pfam fallback for ko-less targets (e.g. archaeal amoA / PF12942).
        for pf in pfams:
            tmp = PFAM_PINNED / f"{pf}.hmm"
            if tmp.exists() and not args.upstream:
                print(f"  + pinned Pfam fallback: {pf}")
            else:
                (CACHE / "ko_hmms").mkdir(parents=True, exist_ok=True)
                tmp = CACHE / "ko_hmms" / f"{pf}.hmm"
                if not fetch_hmm_interpro(pf, tmp):
                    continue
                print(f"  + fetched Pfam fallback: {pf}")
            out.write(tmp.read_text().rstrip("\n") + "\n")
            tc_rows.append(f"{pf}\t\t{extract_tc(tmp) or ''}")

    TC_TSV.write_text("\n".join(tc_rows) + "\n")
    import subprocess
    subprocess.run(["hmmpress", "-f", str(CONCAT)], check=True)
    print(f"[build_hmm_db] wrote {CONCAT} ({len(found)} KO + custom + Pfam profiles) and pressed indices")


if __name__ == "__main__":
    main()
