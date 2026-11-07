#!/usr/bin/env bash
# Launcher input checks (review 2026-11, phase R1): each case must be refused before
# Snakemake starts, with a message that names the file; the smoke panel must still plan
# N samples. Dry runs only, no pipeline job. Run inside the conda env:
#   validation/launcher_cases/run_cases.sh [panel_dir] [n_samples]      (default: test_data 5)
set -u
cd "$(dirname "$0")/../.."
PANEL=${1:-test_data}; NSAMPLES=${2:-5}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
fail=0
run() { python scycle.py --input "$1" --dry-run --plain --skip-db-setup --prodigal-mode single 2>&1; }
expect() {  # expect <case> <exit> <pattern>: run on $T/<case>, check exit status and output
  out=$(run "$T/$1"); rc=$?
  if [ "$rc" -eq "$2" ] && grep -q -- "$3" <<<"$out"; then echo "ok    $1"
  else echo "FAIL  $1 (exit $rc, wanted $2 and /$3/)"; echo "$out" | tail -5; fail=1; fi
}
mkdir -p "$T/space" "$T/dup" "$T/notfasta" "$T/badchar"
printf '>p1\nMKVLAAGIVALLLAAGCSS\n' > "$T/space/bad name.faa"
printf '>p1\nMKVLAAGIVALLLAAGCSS\n' > "$T/dup/a.faa"; printf '>c1\nATGAAAGTTCTGGCA\n' > "$T/dup/a.fna"
printf 'not a fasta\nMKVLAAGIV\n' > "$T/notfasta/x.faa"
printf '>p1\nMKVL?AGIVALLL\n' > "$T/badchar/q.faa"
expect space    1 "bad name"
expect dup      1 "a.faa and .*a.fna"
expect notfasta 1 "not a FASTA file"
expect badchar  1 "illegal character '?'"
out=$(run "$PANEL"); rc=$?
if [ "$rc" -eq 0 ] && grep -q "$NSAMPLES sample(s)" <<<"$out"; then echo "ok    $PANEL: dry run plans $NSAMPLES samples"
else echo "FAIL  $PANEL (exit $rc)"; echo "$out" | tail -5; fail=1; fi
exit $fail
