#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

python - "$ROOT" "$TMPDIR" <<'PY'
import os
import re
import sys
from pathlib import Path

from propka.run import single


ROOT = Path(sys.argv[1])
TMPDIR = Path(sys.argv[2])
PDB_DIR = ROOT / "tests" / "pdb"
RESULTS_DIR = ROOT / "tests" / "results"


def parse_pka(path):
    values = []
    at_pka = False
    with open(path, "rt") as handle:
        for line in handle:
            if at_pka:
                if line.startswith("---"):
                    at_pka = False
                else:
                    match = re.search(r"\d+\.\d+", line[13:])
                    if match is None:
                        raise ValueError(
                            "Could not parse pKa value from {0:s}: {1:s}".format(
                                str(path), line.rstrip()))
                    values.append(float(match.group()))
            elif "model-pKa" in line:
                at_pka = True
    return values


def write_dat(path, values):
    with open(path, "wt") as handle:
        for value in values:
            handle.write("{0:8.2f}\n".format(value))


def read_dat(path):
    return [float(line) for line in path.read_text().splitlines()]


def run_case(stem):
    input_path = PDB_DIR / "{0:s}.pdb".format(stem)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    workdir = TMPDIR / stem
    workdir.mkdir()
    cwd = Path.cwd()
    try:
        os.chdir(workdir)
        single(str(input_path.resolve()))
        return parse_pka(workdir / "{0:s}.pka".format(stem))
    finally:
        os.chdir(cwd)


def mismatches(left, right, tolerance=0.01):
    return [
        (index, a, b, a - b)
        for index, (a, b) in enumerate(zip(left, right))
        if abs(a - b) > tolerance
    ]


outputs = [
    ("3SGB", "3SGB_new.dat"),
    ("3SGB_noicode", "3SGB_new_noicode.dat"),
]

for stem, result_name in outputs:
    values = run_case(stem)
    result_path = RESULTS_DIR / result_name
    write_dat(result_path, values)
    print("{0:s}: wrote {1:d} pKa values to {2:s}".format(
        stem, len(values), str(result_path.relative_to(ROOT))))

old_values = read_dat(RESULTS_DIR / "3SGB.dat")
new_values = read_dat(RESULTS_DIR / "3SGB_new.dat")
noicode_values = read_dat(RESULTS_DIR / "3SGB_new_noicode.dat")

if len(new_values) != len(old_values):
    raise SystemExit("3SGB_new.dat length differs from 3SGB.dat")
if len(noicode_values) != len(old_values):
    raise SystemExit("3SGB_new_noicode.dat length differs from 3SGB.dat")

noicode_mismatch = mismatches(noicode_values, old_values)
if noicode_mismatch:
    print("3SGB_new_noicode.dat does not match legacy 3SGB.dat:")
    for item in noicode_mismatch:
        print("  index {0:d}: new_noicode={1:.2f}, legacy={2:.2f}, diff={3:.2f}".format(
            item[0], item[1], item[2], item[3]))
    raise SystemExit(1)

new_mismatch = mismatches(new_values, old_values)
print("3SGB_new_noicode.dat matches legacy 3SGB.dat")
print("3SGB_new.dat differs from legacy 3SGB.dat at {0:d} pKa values".format(
    len(new_mismatch)))
for item in new_mismatch:
    print("  index {0:d}: new={1:.2f}, legacy={2:.2f}, diff={3:.2f}".format(
        item[0], item[1], item[2], item[3]))
PY
