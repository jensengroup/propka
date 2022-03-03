import propka.atom
import propka.protonate


def test_protonate_atom():
    atom = propka.atom.Atom(
        "HETATM 4479  V   VO4 A1578     -19.097  16.967   0.500  1.00 17.21           V  "
    )
    assert not atom.is_protonated
    p = propka.protonate.Protonate()
    p.protonate_atom(atom)
    assert atom.is_protonated
    assert atom.number_of_protons_to_add == 6
