import propka.atom
import propka.conformation_container
import propka.hydrogens
import propka.parameters
import propka.protonate
import propka.vector_algebra


def test_protonate_atom():
    atom = propka.atom.Atom(
        "HETATM 4479  V   VO4 A1578     -19.097  16.967   0.500  1.00 17.21           V  "
    )
    assert not atom.is_protonated
    p = propka.protonate.Protonate()
    p.protonate_atom(atom)
    assert atom.is_protonated
    assert atom.number_of_protons_to_add == 6


def make_atom_with_insertion_code():
    """Create an atom attached to a container for protonation tests."""
    atom = propka.atom.Atom()
    atom.set_property(
        name="CA", res_name="GLY", chain_id="XX", res_num=10, icode="A")
    atom.element = "C"
    atom.type = "atom"
    atom.number_of_protons_to_add = 1
    container = propka.conformation_container.ConformationContainer(
        "test", propka.parameters.Parameters(), None)
    container.add_atom(atom)
    return atom


def test_protonate_preserves_insertion_code_on_generated_hydrogen():
    """The general protonator copies insertion codes to new hydrogens."""
    atom = make_atom_with_insertion_code()

    propka.protonate.Protonate.add_proton(
        atom, propka.vector_algebra.Vector(1.0, 2.0, 3.0))

    hydrogen = atom.bonded_atoms[0]
    assert hydrogen.residue_key == atom.residue_key
    assert hydrogen.residue_label.endswith("XX:A")


def test_legacy_protonate_preserves_insertion_code_on_hydrogen():
    """The legacy protonator also copies insertion codes to hydrogens."""
    atom = make_atom_with_insertion_code()

    hydrogen = propka.hydrogens.make_new_h(atom, 1.0, 2.0, 3.0)

    assert hydrogen.residue_key == atom.residue_key
    assert hydrogen.residue_label.endswith("XX:A")
