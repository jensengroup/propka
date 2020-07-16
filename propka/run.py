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


def single(filename: str, optargs: list, stream=None, write_pka: bool = True):
    """Run a single PROPKA calculation using ``filename`` as input.

    Args:
        filename (str): name of input file. If filestream is not passed via
            ``stream``, should be path to the file to be read.
        optargs (list): Optional, commandline options for propka.
        stream : optional filestream handle. If ``None``, then ``filename``
            will be used as path to input file for reading.
        write_pka (bool): Controls if the pKa file should be writen to disk.

    Example:
    Given an input file "protein.pdb", run the equivalent of ``propka3
    --mutation=N25R/N181D -v --pH=7.2 protein.pdb`` as::

        propka.run.single("protein.pdb",
                          optargs=["--mutation=N25R/N181D", "-v", "--pH=7.2"])

    By default, a pKa file will be written. However in some cases one may wish
    to not output this file and just have access to the
    MolecularContainer object. If so, then pass ``False`` to ``write_pka``::

        mol = propka.run.single("protein.pdb", write_pka=False)

    In some cases, one may also want to pass a file-like (e.g. StringIO) object
    instead of a file path as a string. In these cases the file-like object
    should be passed to the ``stream`` argument. Since file-like objects do not
    usually have names, and the ``filename`` argument must be provided to
    help propka determine the input file type, and assigns the file name for
    the MolecularContainer object::

        mol = propka.run.single('input.pdb', stream=string_io_file)

    In this case, a PDB file-like object was passed as `string_io_file`. The
    resultant pKa file will be written out as `input.pka`.

    """
    # Deal with input optarg options
    optargs = optargs if optargs is not None else []
    optargs += [filename]
    options = loadOptions(optargs)

    parameters = read_parameter_file(options.parameters, Parameters())

    # Only filename present should be the one passed via the arguments list
    # Anything else will probably have been passed using optargs and the `-f`
    # flag. Note: can be mulitple -f entries and all are appended before the
    # input filename
    if len(options.filenames) > 1:
        ignored_list = options.filenames[:-1]
        # overwrite the options.filenames entry
        options.filenames = options.filenames[-1]
        _LOGGER.warning(f"Ignoring filenames: {ignored_list}")

    my_molecule = MolecularContainer(parameters, options)
    my_molecule = read_molecule_file(filename, my_molecule, stream=stream)
    my_molecule.calculate_pka()

    # write outputs
    if options.generate_propka_input:
        my_molecule.write_propka()
    if write_pka:
        my_molecule.write_pka()

    return my_molecule
