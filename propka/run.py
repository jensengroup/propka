"""Entry point for PROPKA script."""
import logging
from propka.lib import loadOptions
from propka.molecular_container import Molecular_container


_LOGGER = logging.getLogger("PROPKA")


def main(optargs=None):
    """Read in structure files, calculate pKa values, and print pKa files."""
    # loading options, flags and arguments
    optargs = optargs if optargs is not None else []
    options = loadOptions(*optargs)
    pdbfiles = options.filenames
    for pdbfile in pdbfiles:
        my_molecule = Molecular_container(pdbfile, options)
        my_molecule.calculate_pka()
        my_molecule.write_pka()


def single(pdbfile, optargs=None):
    """Run a single PROPKA calculation using *pdbfile* as input.

    Commandline options can be passed as a **list** in *optargs*.

    .. rubric:: Example

    ::
       single("protein.pdb", optargs=["--mutation=N25R/N181D", "-v", "--pH=7.2"])
    """
    optargs = optargs if optargs is not None else []
    options = loadOptions(*optargs)
    pdbfile = options.filenames.pop(0)
    if len(options.filenames) > 0:
        _LOGGER.warning("Ignoring filenames: %s", options.filenames)
    my_molecule = Molecular_container(pdbfile, options)
    my_molecule.calculate_pka()
    my_molecule.write_pka()
    return my_molecule
