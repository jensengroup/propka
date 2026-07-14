import zipfile
from io import StringIO

import pytest

import propka.input as m
from propka.lib import loadOptions
from propka.molecular_container import MolecularContainer
from propka.parameters import Parameters


def test_open_file_for_reading(tmp_path):
    path = tmp_path / "tmp.txt"
    path.write_text("One\nTwo\nThree\n")
    # str
    with m.open_file_for_reading(str(path)) as outer:
        assert outer.read() == "One\nTwo\nThree\n"
    assert outer.closed
    # Path
    with m.open_file_for_reading(path) as outer:
        # TextIO
        with m.open_file_for_reading(outer) as inner:
            assert inner.readline() == "One\n"
        assert not outer.closed
        assert outer.readline() == "Two\n"
    assert outer.closed


def test_open_file_for_reading__zipfile(tmp_path):
    zippath = tmp_path / "tmp.zip"
    arcname = "foo/bar.txt"
    with zipfile.ZipFile(zippath, "w") as ziphandle:
        ziphandle.writestr(arcname, "One\nTwo\nThree\n")
    with m.open_file_for_reading(zippath / arcname) as outer:
        assert outer.readline() == "One\n"
    assert outer.closed


def test_read_mmcif_atom_site_identity_fields():
    cif = StringIO("""data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 100001 N N . GLY L1 10000 B 1.0 2.0 3.0 1.0 10.0 10000 GLY AA N 1
ATOM 100002 C CA . GLY L1 10000 B 2.0 3.0 4.0 1.0 11.0 10000 GLY AA CA 1
""")
    options = loadOptions(["input.cif"])
    parameters = Parameters()
    molecule = MolecularContainer(parameters, options)

    conformations, names = m.read_mmcif(cif, parameters, molecule)

    assert names == ["1A"]
    atoms = conformations["1A"].atoms
    assert [atom.numb for atom in atoms] == [100001, 100002]
    assert {atom.chain_id for atom in atoms} == {"AA"}
    assert {atom.residue_key for atom in atoms} == {("AA", 10000, "B")}
    assert atoms[0].terminal == "N+"
    assert atoms[0].element == "N"


def test_read_mmcif_missing_required_atom_site_value():
    cif = StringIO("""data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . GLY A 1 ? 1.0 2.0 1.0 10.0
""")
    options = loadOptions(["input.cif"])
    parameters = Parameters()
    molecule = MolecularContainer(parameters, options)

    with pytest.raises(ValueError, match="Cartn_z"):
        m.read_mmcif(cif, parameters, molecule)


def test_read_mmcif_invalid_optional_model_number():
    cif = StringIO("""data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 ? 1.0 2.0 3.0 1.0 10.0 bad
""")
    options = loadOptions(["input.cif"])
    parameters = Parameters()
    molecule = MolecularContainer(parameters, options)

    with pytest.raises(ValueError, match="pdbx_PDB_model_num"):
        m.read_mmcif(cif, parameters, molecule)


@pytest.mark.parametrize("altloc", ["", "AB", "10", "*"])
def test_normalize_altloc_rejects_unsupported_identifiers(altloc):
    with pytest.raises(ValueError, match="alternate-location"):
        m._normalize_altloc(altloc)
