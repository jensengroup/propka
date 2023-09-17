from propka.input import read_parameter_file, read_molecule_file
from propka.lib import loadOptions
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
