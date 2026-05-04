#!/usr/bin/env python3
import argparse
import shutil
import subprocess
from pathlib import Path


def run(cmd):
    return subprocess.run(cmd, check=False, capture_output=True, text=True)


def get_id(fmt, path):
    # Prefer InChI (more canonical), fall back to SMILES
    res = run(["obabel", f"-i{fmt}", str(path), "-oinchi"])
    if res.returncode == 0:
        for line in res.stdout.splitlines():
            line = line.strip()
            if line.startswith("InChI="):
                return line

    res = run(["obabel", f"-i{fmt}", str(path), "-osmi"])
    if res.returncode == 0:
        for line in res.stdout.splitlines():
            line = line.strip()
            if line:
                return line.split()[0]
    return ""


def main():
    parser = argparse.ArgumentParser(
        description="Convert Gaussian16 .log files to PDB and verify connectivity."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Root directory to search (default: current directory).",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output directory (default: <root>/pdb_from_g16).",
    )
    args = parser.parse_args()

    if shutil.which("obabel") is None:
        raise SystemExit("ERROR: obabel not found. Install Open Babel (obabel).")

    root = Path(args.root).resolve()
    out_dir = Path(args.out).resolve() if args.out else root / "pdb_from_g16"
    out_dir.mkdir(parents=True, exist_ok=True)

    logs = sorted(root.rglob("*.log"))
    if not logs:
        raise SystemExit(f"No .log files found under: {root}")

    report = out_dir / "verification_report.tsv"
    with report.open("w", encoding="utf-8") as f:
        f.write("log\tpdb\tmatch\tlog_id\tpdb_id\n")
        for log in logs:
            pdb = out_dir / f"{log.stem}.pdb"

            res = run(["obabel", "-ig16", str(log), "-opdb", "-O", str(pdb)])
            if res.returncode != 0:
                f.write(f"{log}\t{pdb}\tfailed\t\t\n")
                continue

            log_id = get_id("g16", log)
            pdb_id = get_id("pdb", pdb)
            match = "ok" if (log_id and pdb_id and log_id == pdb_id) else "mismatch"
            f.write(f"{log}\t{pdb}\t{match}\t{log_id}\t{pdb_id}\n")

    print(f"Done. PDBs written to: {out_dir}")
    print(f"Verification report: {report}")


if __name__ == "__main__":
    main()
