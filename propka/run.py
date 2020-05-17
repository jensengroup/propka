# entry point for propka script
import logging
import propka.lib, propka.molecular_container


_LOGGER = logging.getLogger("PROPKA")


def main():
    """
    Reads in structure files, calculates pKa values, and prints pKa files
    """
    # loading options, flaggs and arguments
    options = propka.lib.loadOptions()
    pdbfiles = options.filenames

    for pdbfile in pdbfiles:
        my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
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
    options = propka.lib.loadOptions(*optargs)
    pdbfile = options.filenames.pop(0)
    if len(options.filenames) > 0:
        _LOGGER.warning("Ignoring filenames: %s", options.filenames)

    my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
    my_molecule.calculate_pka()
    my_molecule.write_pka()
    return my_molecule
