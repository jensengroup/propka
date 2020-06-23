"""
Script functionality
====================

The :mod:`run` module provides a high-level interface to PROPKA 3.

The :program:`propka3` script consists of the :func:`main`
function. If similar functionality is desired from a Python script
(without having to call the :program:`propka` script itself) then the
:func:`single` function can be used instead.

"""
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

    Example
    -------
    Given an input file "protein.pdb", run the equivalent of ``propka3
    --mutation=N25R/N181D -v --pH=7.2 protein.pdb`` as::

       propka.run.single("protein.pdb",
                         optargs=["--mutation=N25R/N181D", "-v", "--pH=7.2"])


    .. todo::
       Test :func:`single`, not sure if it is correctly processing ``pdbfile``.

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
