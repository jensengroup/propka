"""Tests for PROPKA"""
import logging
import os
import re
import json
from pathlib import Path
import pytest
from pytest import approx
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer
from propka.input import read_parameter_file, read_molecule_file
from propka.lib import loadOptions
from typing import List


_LOGGER = logging.getLogger(__name__)


# Number of decimal places for maximum tolerable error.  Set by number of
# decimal places in pKa output as well as need to make unmodified code work
# on WSL Ubuntu 18.04
MAX_ERR_DECIMALS = 2
MAX_ERR_ABS = 10**-MAX_ERR_DECIMALS


# This directory
TEST_DIR = Path("tests")
# Location for test PDBs
PDB_DIR = Path("pdb")
# Location for results for comparing output (allow running from tests/ and
# ../tests/)
RESULTS_DIR = Path("tests/results")
if not RESULTS_DIR.is_dir():
    _LOGGER.warning("Switching to sub-directory")
    RESULTS_DIR = Path("results")


def get_test_dirs():
    """Get locations of test files.

    Returns:
        dictionary with test file locations.
    """
    path_dict = {}
    for key, path in [("pdbs", PDB_DIR), ("results", RESULTS_DIR)]:
        test_path = TEST_DIR / path
        if test_path.is_dir():
            path_dict[key] = test_path
        else:
            test_path = path
            if test_path.is_dir():
                path_dict[key] = test_path
            else:
                errstr = (
                    "Can't find {0:s} test files in {1:s}".format(
                        key, [TEST_DIR / path, path]))
                raise FileNotFoundError(errstr)
    return path_dict


def run_propka(options, pdb_path, tmp_path):
    """Run PROPKA software.

    Args:
        options:  list of PROPKA options
        pdb_path:  path to PDB file
        tmp_path:  path for working directory
    """
    options += [str(pdb_path)]
    args = loadOptions(options)
    cwd = Path.cwd()
    try:
        _LOGGER.warning(
            "Working in tmpdir {0:s} because of PROPKA file output; "
            "need to fix this.".format(str(tmp_path)))
        os.chdir(tmp_path)
        parameters = read_parameter_file(args.parameters, Parameters())
        molecule = MolecularContainer(parameters, args)
        molecule = read_molecule_file(str(pdb_path), molecule)
        molecule.calculate_pka()
        molecule.write_pka()
    finally:
        os.chdir(cwd)


def parse_pka(pka_path: Path) -> dict:
    """Parse testable data from a .pka file into a dictionary.
    """
    pka_list: List[float] = []
    data: dict = {"pKa": pka_list}

    with open(pka_path, "rt") as pka_file:
        at_pka = False
        for line in pka_file:
            if at_pka:
                if line.startswith("---"):
                    at_pka = False
                else:
                    m = re.search(r'\d+\.\d+', line[13:])
                    assert m is not None
                    pka_list.append(float(m.group()))
            elif "model-pKa" in line:
                at_pka = True
            else:
                m = re.match(
                    r"The pI is *(\d+\.\d+) .folded. and *(\d+\.\d+) .unfolded.",
                    line)
                if m is not None:
                    data["pI_folded"] = float(m.group(1))
                    data["pI_unfolded"] = float(m.group(2))

    return data


def compare_output(pdb, tmp_path, ref_path):
    """Compare results of test with reference.

    Args:
        pdb:  PDB filename stem
        tmp_path:  temporary directory
        ref_path:  path with reference results
    Raises:
        ValueError if results disagree.
    """
    with open(ref_path, "rt") as ref_file:
        if ref_path.name.endswith(".json"):
            ref_data = json.load(ref_file)
        else:
            ref_data = {"pKa": [float(line) for line in ref_file]}

    test_data = parse_pka(tmp_path / f"{pdb}.pka")

    for key in ref_data:
        assert test_data[key] == approx(ref_data[key], abs=MAX_ERR_ABS), key


@pytest.mark.parametrize("pdb, options", [
    pytest.param('sample-issue-140', [], id="sample-issue-140: no options"),
    pytest.param("1FTJ-Chain-A", [], id="1FTJ-Chain-A: no options"),
    pytest.param('1HPX', [], id="1HPX: no options"),
    pytest.param('4DFR', [], id="4DFR: no options"),
    pytest.param('3SGB', [], id="3SGB: no options"),
    pytest.param('3SGB-subset', [
        "--titrate_only",
        "E:17,E:18,E:19,E:29,E:44,E:45,E:46,E:118,E:119,E:120,E:139"],
                 id="3SGB: --titrate_only"),
    pytest.param('1HPX-warn', ['--quiet'], id="1HPX-warn: --quiet")])
def test_regression(pdb, options, tmp_path):
    """Basic regression test of PROPKA functionality."""
    path_dict = get_test_dirs()
    ref_path = None

    for ext in ["json", "dat"]:
        ref_path = path_dict["results"] / f"{pdb}.{ext}"
        if ref_path.is_file():
            ref_path = ref_path.resolve()
            break
    else:
        _LOGGER.warning("Missing results file for comparison: {0:s}".format(
            str(ref_path)))
        ref_path = None
    pdb_path = path_dict["pdbs"] / ("{0:s}.pdb".format(pdb))
    if pdb_path.is_file():
        pdb_path = pdb_path.resolve()
    else:
        errstr = "Missing PDB file: {0:s}".format(pdb_path)
        raise FileNotFoundError(errstr)
    tmp_path = Path(tmp_path).resolve()

    run_propka(options, pdb_path, tmp_path)
    if ref_path is not None:
        compare_output(pdb, tmp_path, ref_path)
