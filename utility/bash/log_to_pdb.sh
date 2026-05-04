#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar

root="${1:-.}"
out_dir="${2:-$root/pdb_from_g16}"

if ! command -v obabel >/dev/null 2>&1; then
  echo "ERROR: obabel not found. Install Open Babel (obabel)." >&2
  exit 1
fi

mkdir -p "$out_dir"

logs=( "$root"/**/*.log )
if (( ${#logs[@]} == 0 )); then
  echo "No .log files found under: $root" >&2
  exit 1
fi

report="$out_dir/verification_report.tsv"
printf "log\tpdb\tmatch\tlog_id\tpdb_id\n" > "$report"

get_id() {
  local fmt="$1"
  local file="$2"
  local out=""
  local id=""

  if out=$(obabel -i"$fmt" "$file" -oinchi 2>/dev/null); then
    id=$(printf "%s\n" "$out" | awk '/^InChI=/{print; exit}')
    if [[ -n "$id" ]]; then
      printf "%s\n" "$id"
      return 0
    fi
  fi

  out=$(obabel -i"$fmt" "$file" -osmi 2>/dev/null || true)
  id=$(printf "%s\n" "$out" | awk 'NF{print $1; exit}')
  if [[ -n "$id" ]]; then
    printf "%s\n" "$id"
    return 0
  fi

  return 1
}

for log in "${logs[@]}"; do
  base="$(basename "$log" .log)"
  pdb="$out_dir/$base.pdb"

  if ! obabel -ig16 "$log" -opdb -O "$pdb" >/dev/null 2>&1; then
    echo "Failed to convert: $log" >&2
    continue
  fi

  log_id="$(get_id g16 "$log" || true)"
  pdb_id="$(get_id pdb "$pdb" || true)"

  if [[ -n "$log_id" && -n "$pdb_id" && "$log_id" == "$pdb_id" ]]; then
    match="ok"
  else
    match="mismatch"
  fi

  printf "%s\t%s\t%s\t%s\t%s\n" "$log" "$pdb" "$match" "$log_id" "$pdb_id" >> "$report"
done

echo "Done. PDBs written to: $out_dir"
echo "Verification report: $report"
