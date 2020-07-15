"""Tests for PROPKA stream io"""
import logging
from pathlib import Path
from io import StringIO
import pytest
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer
from propka.input import read_parameter_file, read_molecule_file
from propka.lib import loadOptions

from .test_basic_regression import get_test_dirs, compare_output


_LOGGER = logging.getLogger(__name__)


def get_paths(pdb):
    """Helper function to get the path to the input and reference files"""
    path_dict = get_test_dirs()
    ref_path = path_dict["results"] / ("{0:s}.dat".format(pdb))
    pdb_path = path_dict["pdbs"] / ("{0:s}.pdb".format(pdb))

    return ref_path.resolve(), pdb_path.resolve()


def run_propka_stream(options, input_file, filename):
    """Run PROPKA software.

    Args:
        options:  list of PROPKA options
        input_file:  file-like PDB object
        filename: filename for the file-like PDB object
        tmp_path:  path for working directory
    """
    options += [filename]
    args = loadOptions(options)
    parameters = read_parameter_file(args.parameters, Parameters())
    molecule = MolecularContainer(parameters, args)
    molecule = read_molecule_file(input_file, molecule, filename)
    molecule.calculate_pka()
    molecule.write_pka()
    if args.generate_propka_input:
        molecule.write_propka()


@pytest.mark.parametrize("pdb, options", [
    pytest.param("1FTJ-Chain-A", [], id="1FTJ-Chain-A: no options"),
    pytest.param('3SGB-subset', [
        "--titrate_only",
        "E:17,E:18,E:19,E:29,E:44,E:45,E:46,E:118,E:119,E:120,E:139"],
                 id="3SGB: --titrate_only"),
    pytest.param('1HPX-warn', ['--quiet'], id="1HPX-warn: --quiet"),
])
def test_textio_filestream(tmpdir, pdb, options):
    """Basic regression test using TextIO streams for the input PDB file"""
    # Get the relevant paths
    ref_path, pdb_path = get_paths(pdb)
    filename = f"{pdb}.pdb"

    filestream = open(pdb_path, 'r')

    with tmpdir.as_cwd():
        run_propka_stream(options, filestream, filename)
        compare_output(pdb, Path.cwd(), ref_path)

    filestream.close()


@pytest.mark.parametrize("pdb, options", [
    pytest.param("1FTJ-Chain-A", [], id="1FTJ-Chain-A: no options"),
    pytest.param('3SGB-subset', [
        "--titrate_only",
        "E:17,E:18,E:19,E:29,E:44,E:45,E:46,E:118,E:119,E:120,E:139"],
                 id="3SGB: --titrate_only"),
    pytest.param('1HPX-warn', ['--quiet'], id="1HPX-warn: --quiet"),
])
def test_stringio_filestream(tmpdir, pdb, options):
    """Basic regression test using StringIO streams for the input PDB file"""
    # Get the relevant paths
    ref_path, pdb_path = get_paths(pdb)
    filename = f"{pdb}.pdb"

    with open(pdb_path, 'r') as writer:
        filestream = StringIO(writer.read())

    with tmpdir.as_cwd():
        run_propka_stream(options, filestream, filename)
        compare_output(pdb, Path.cwd(), ref_path)

    filestream.close()


def test_typerror_nofilename(tmpdir):
    """Tests for raised TypeError when not passing a filename to
    read_molecule_file and using a file-like object without a name"""
    pdb = "1FTJ-Chain-A"
    options = []

    ref_path, pdb_path = get_paths(pdb)

    with open(pdb_path, 'r') as writer:
        filestream = StringIO(writer.read())

    with tmpdir.as_cwd():
        errmsg = "Path of provided input_file could not be determined"
        with pytest.raises(TypeError, match=errmsg):
            # default value of filename is None
            run_propka_stream(options, filestream, filename=None)
