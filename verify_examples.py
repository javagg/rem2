#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
WAVEPORT_CONFIG = REPO_ROOT / "examples" / "cpw" / "waveport_tri3_verify.json"
WAVEPORT_OUTPUT = REPO_ROOT / "output" / "waveport_tri3_verify_smoke"


def run_waveport_smoke() -> int:
    if not WAVEPORT_CONFIG.exists():
        print(f"missing config: {WAVEPORT_CONFIG}", file=sys.stderr)
        return 1

    if WAVEPORT_OUTPUT.exists():
        shutil.rmtree(WAVEPORT_OUTPUT)

    cmd = [
        "cargo",
        "run",
        "-p",
        "rem-cli",
        "--",
        str(WAVEPORT_CONFIG.relative_to(REPO_ROOT)),
        "-o",
        str(WAVEPORT_OUTPUT.relative_to(REPO_ROOT)),
    ]
    print("running:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        return result.returncode

    postpro = WAVEPORT_OUTPUT / "postpro"
    support_csv = postpro / "wave-port-support.csv"
    peak_csv = postpro / "domain-E-peak-by-tag.csv"

    missing = [path for path in (support_csv, peak_csv) if not path.exists()]
    if missing:
        for path in missing:
            print(f"missing output: {path}", file=sys.stderr)
        return 2

    support_rows = read_csv_rows(support_csv)
    peak_rows = read_csv_rows(peak_csv)
    if not support_rows:
        print(f"no support-region rows in {support_csv}", file=sys.stderr)
        return 3
    if not peak_rows:
        print(f"no peak-domain-energy rows in {peak_csv}", file=sys.stderr)
        return 4

    positive_energies = []
    for row in peak_rows:
        try:
            energy = float(row["Electric Field Energy (J)"])
        except (KeyError, ValueError) as exc:
            print(f"bad peak energy row {row}: {exc}", file=sys.stderr)
            return 5
        if energy > 0.0:
            positive_energies.append(energy)

    if not positive_energies:
        print(f"all peak-domain energies are zero in {peak_csv}", file=sys.stderr)
        return 6

    print(f"verified {support_csv.relative_to(REPO_ROOT)} with {len(support_rows)} row(s)")
    print(
        f"verified {peak_csv.relative_to(REPO_ROOT)} with {len(peak_rows)} row(s); "
        f"max energy = {max(positive_energies):.6e} J"
    )
    return 0


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="REM example verification helpers")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser(
        "waveport-smoke",
        help="run the minimal CPW WavePort example and validate driven diagnostics CSVs",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "waveport-smoke":
        return run_waveport_smoke()
    parser.error(f"unsupported command: {args.command}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())