"""Tests for PROPKA's run module"""
import logging
import os
from pathlib import Path
from io import StringIO
import pytest
import propka.run as pkrun

from .test_basic_regression import compare_output
from .test_streamio import get_paths


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb, options", [
    pytest.param("1FTJ-Chain-A", (), id="1FTJ-Chain-A: no options"),
    pytest.param('3SGB-subset', (
        "--titrate_only",
        "E:17,E:18,E:19,E:29,E:44,E:45,E:46,E:118,E:119,E:120,E:139"),
                 id="3SGB: --titrate_only"),
    pytest.param('1HPX-warn', ('--quiet',), id="1HPX-warn: --quiet"),
])
def test_single_file(tmpdir, pdb, options):
    """Basic regression test using propka.run.single and local file for the
    input PDB file"""
    ref_path, pdb_path = get_paths(pdb)
    filename = str(pdb_path)

    with tmpdir.as_cwd():
        pkrun.single(filename, options)
        compare_output(pdb, Path.cwd(), ref_path)
        assert os.path.isfile(f'{pdb}.pka')


@pytest.mark.parametrize("pdb, options", [
    pytest.param("1FTJ-Chain-A", (), id="1FTJ-Chain-A: no options"),
    pytest.param('3SGB-subset', (
        "--titrate_only",
        "E:17,E:18,E:19,E:29,E:44,E:45,E:46,E:118,E:119,E:120,E:139"),
                 id="3SGB: --titrate_only"),
    pytest.param('1HPX-warn',('--quiet',), id="1HPX-warn: --quiet"),
])
def test_single_filestream(tmpdir, pdb, options):
    """Basic regression test using StringIO streams for the input PDB file"""
    ref_path, pdb_path = get_paths(pdb)
    filename = f"{pdb}.pdb"

    with open(pdb_path, 'r') as writer:
        filestream = StringIO(writer.read())

    with tmpdir.as_cwd():
        pkrun.single(filename, options, stream=filestream)
        compare_output(pdb, Path.cwd(), ref_path)
        assert os.path.isfile(f'{pdb}.pka')

    filestream.close()


def test_single_nopka(tmpdir):
    """Basic test to check that the pKa file is not written when write_pka is
    `False`"""
    pdb = "1FTJ-Chain-A"
    ref_path, pdb_path = get_paths(pdb)
    filename = f"{pdb}.pdb"

    with open(pdb_path, 'r') as writer:
        filestream = StringIO(writer.read())

    pkrun.single(filename, stream=filestream, write_pka=False)
    assert not os.path.isfile(f"{pdb}.pka")


def test_single_extra_files_logwarn(tmpdir, caplog):
    """Tests that a logging warning is thrown if passing files via optargs"""
    pdb = "1FTJ-Chain-A"
    options = ('-f foo.pdb bar.pdb', '-f test.pdb test2.pdb')
    ref_path, pdb_path = get_paths(pdb)
    filename = str(pdb_path)

    with tmpdir.as_cwd():
        pkrun.single(filename, options)

        wmsg = ("Ignoring extra filenames passed: [' foo.pdb bar.pdb', "
                "' test.pdb test2.pdb']")
        assert wmsg in caplog.records[0].message
