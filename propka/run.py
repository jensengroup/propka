"""Entry point for PROPKA script."""
import logging
from propka.lib import loadOptions
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer


_LOGGER = logging.getLogger("PROPKA")


def main(optargs=None):
    """Read in structure files, calculate pKa values, and print pKa files."""
    # loading options, flags and arguments
    optargs = optargs if optargs is not None else []
    options = loadOptions(*optargs)
    pdbfiles = options.filenames
    parameters = read_parameter_file(options.parameters, Parameters())
    for pdbfile in pdbfiles:
        my_molecule = MolecularContainer(parameters, options)
        my_molecule = read_molecule_file(pdbfile, my_molecule)
        my_molecule.calculate_pka()
        my_molecule.write_pka()
        if options.generate_propka_input:
            my_molecule.write_propka()


def single(pdbfile, optargs=None):
    """Run a single PROPKA calculation using *pdbfile* as input.

    Commandline options can be passed as a **list** in *optargs*.

    .. rubric:: Example

    ::
       single("protein.pdb", optargs=["--mutation=N25R/N181D", "-v",
              "--pH=7.2"])
    """
    optargs = optargs if optargs is not None else []
    options = loadOptions(*optargs)
    pdbfile = options.filenames.pop(0)
    parameters = read_parameter_file(options.parameters, Parameters())
    if len(options.filenames) > 0:
        _LOGGER.warning("Ignoring filenames: {0:s}".format(options.filenames))
    my_molecule = MolecularContainer(parameters, options)
    my_molecule = read_molecule_file(pdbfile, my_molecule)
    my_molecule.calculate_pka()
    my_molecule.write_pka()
    if options.generate_propka_input:
        my_molecule.write_propka()
    return my_molecule
