#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar
root="${1:-.}"
out_dir="${2:-$root/smi_from_g16}" # default dir
trim="${3:-$root}"
trim="${trim%/}"
mkdir -p "$out_dir"
for f in "$root"/**/*.log; do
  rel="${f#"$trim"/}"
  rel="/${rel#/}"            # -> /qchem_data/log_to_pdb_in/e1.log
  base="$(basename "$f" .log)"
  out="$out_dir/$base.smi"
  obabel -ig16 "$f" -ocan | awk -v t="$rel" 'NR==1{print $1 "\t" t}' > "$out"
done
