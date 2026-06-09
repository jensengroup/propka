from propka.input import read_parameter_file, read_molecule_file
from propka.lib import loadOptions
from propka.atom import Atom
from propka.conformation_container import ConformationContainer
from propka.molecular_container import MolecularContainer
from propka.parameters import Parameters
from pathlib import Path
from pytest import approx

TESTS = Path(__file__).resolve().parent


def load_molecule(code):
    pdb_path = TESTS / f"pdb/{code}.pdb"
    options = [str(pdb_path)]
    args = loadOptions(options)
    parameters = read_parameter_file(args.parameters, Parameters())
    molecule = MolecularContainer(parameters, args)
    molecule = read_molecule_file(str(pdb_path), molecule)
    return molecule


def make_atom(name, chain_id, res_num, icode=" ", res_name="GLY"):
    """
    A small func added to create an object for test
    """
    atom = Atom()
    atom.set_property(
        numb=1,
        name=name,
        res_name=res_name,
        chain_id=chain_id,
        res_num=res_num,
        icode=icode,
        x=0.0,
        y=0.0,
        z=0.0,
        occ="1.0",
        beta="0.0",
    )
    atom.element = name[0]
    atom.type = "atom"
    return atom


def test_ConformationContainer_sort_atoms_multi_character_chain():
    """
    Create multi-chain mmCIF object for test 
    """
    conf = ConformationContainer("test", Parameters(), None)
    atoms = [
        make_atom("CA", "B", 1),
        make_atom("CA", "AA", 1, "B"), #Explict wrong rank to verify sort function
        make_atom("CA", "AA", 1, "A"),
        make_atom("CA", "AA", 2),
    ]
    conf.atoms = atoms

    conf.sort_atoms()

    assert [atom.residue_key for atom in conf.atoms] == [
        ("AA", 1, "A"),
        ("AA", 1, "B"),
        ("AA", 2, " "),
        ("B", 1, " "),
    ]


def test_ConformationContainer_top_up_keeps_insertion_codes_separate():
    """
    To test whether the multi insertion code would be properly dealt with
    """
    conf = ConformationContainer("test", Parameters(), None)
    conf.add_atom(make_atom("CA", "AA", 10, "A", "GLY"))

    ref_atoms = [
        make_atom("CA", "AA", 10, "A", "GLY"),
        make_atom("CA", "AA", 10, "B", "SER"),
    ]
    conf.top_up_from_atoms(ref_atoms)

    assert [(atom.residue_key, atom.res_name) for atom in conf.atoms] == [
        (("AA", 10, "A"), "GLY"),
        (("AA", 10, "B"), "SER"),
    ]


def test_MolecularContainer_get_pi():
    molecule = load_molecule("1HPX")
    molecule.average_of_conformations()
    molecule.calculate_pka()
    pi_folded, pi_unfolded = molecule.get_pi()
    assert pi_folded == approx(9.54, abs=1e-2)
    assert pi_unfolded == approx(8.90, abs=1e-2)


def test_MolecularContainer_top_up_conformations():
    molecule = load_molecule("conf-alt-AB")
    assert len(molecule.conformations) == 2
    assert len(molecule.conformations["1A"]) == 16
    assert len(molecule.conformations["1B"]) == 16

    molecule = load_molecule("conf-alt-BC")
    assert len(molecule.conformations) == 3
    assert len(molecule.conformations["1A"]) == 16
    assert len(molecule.conformations["1B"]) == 16
    assert len(molecule.conformations["1C"]) == 16

    molecule = load_molecule("conf-alt-AB-mutant")
    assert len(molecule.conformations) == 2
    assert len(molecule.conformations["1A"]) == 16
    assert len(molecule.conformations["1B"]) == 17

    molecule = load_molecule("conf-model-mutant")
    assert len(molecule.conformations) == 2
    assert len(molecule.conformations["1A"]) == 16
    assert len(molecule.conformations["2A"]) == 17

    molecule = load_molecule("conf-model-missing-atoms")
    assert len(molecule.conformations) == 3
    assert len(molecule.conformations["1A"]) == 17
    assert len(molecule.conformations["2A"]) == 17
    assert len(molecule.conformations["3A"]) == 17
