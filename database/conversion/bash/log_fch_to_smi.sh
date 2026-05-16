#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar

# Walk $root (recursive **), convert Gaussian logs and formatted checkpoint files
# to one SMILES-per-file in $out_dir (same folder for all types).

root="${1:-.}"
out_dir="${2:-$root/smi_from_g16}"
trim="${3:-$root}"
trim="${trim%/}"

mkdir -p "$out_dir"

to_smi() {
  local f=$1
  local ifmt=$2
  local rel="${f#"$trim"/}"
  rel="/${rel#/}"
  local base
  base="$(basename "${f%.*}")"
  local out="$out_dir/$base.smi"
  obabel -i"$ifmt" "$f" -ocan | awk -v t="$rel" 'NR==1{print $1 "\t" t}' >"$out"
}

# One output directory: every match is handled here (logs + fchk in the same tree).
for f in "$root"/**/*.log "$root"/**/*.fch "$root"/**/*.fchk; do
  case "${f,,}" in
    *.log)  to_smi "$f" "g16" ;;
    *.fch|*.fchk) to_smi "$f" "fch" ;;
  esac
done
