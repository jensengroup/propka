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
import sys
from typing import IO, Iterable, Optional
from propka.lib import loadOptions
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer
from propka.output import _PathArg


_LOGGER = logging.getLogger(__name__)


def main(optargs=None):
    """Read in structure files, calculate pKa values, and print pKa files.


    .. versionchanged:: 3.4.0
       Removed ability to write out PROPKA input files.
    """
    # loading options, flags and arguments
    logger = logging.getLogger("")
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(stdout_handler)
    optargs = optargs if optargs is not None else []
    options = loadOptions(*optargs)
    pdbfiles = options.filenames
    parameters = read_parameter_file(options.parameters, Parameters())
    for pdbfile in pdbfiles:
        my_molecule = MolecularContainer(parameters, options)
        my_molecule = read_molecule_file(pdbfile, my_molecule)
        my_molecule.calculate_pka()
        my_molecule.write_pka()


def single(filename: _PathArg,
           optargs: Iterable[str] = (),
           stream: Optional[IO[str]] = None,
           write_pka: bool = True):
    """Run a single PROPKA calculation using ``filename`` as input.

    Args:
        filename (str): name of input file. If filestream is not passed via
            ``stream``, should be a path to the file to be read.
        optargs (tuple): Optional, commandline options for propka. Extra files
            passed via ``optargs`` will be ignored, see Notes.
        stream : optional filestream handle. If ``None``, then ``filename``
            will be used as path to input file for reading.
        write_pka (bool): Controls if the pKa file should be writen to disk.

    Returns:
        :class:`~propka.molecular_container.MolecularContainer` object.

    Examples:
        Given an input file "protein.pdb", run the equivalent of ``propka3
        --mutation=N25R/N181D -v --pH=7.2 protein.pdb`` as::

            propka.run.single("protein.pdb",
                optargs=["--mutation=N25R/N181D", "-v", "--pH=7.2"])

        By default, a pKa file will be written. However in some cases one may
        wish to not output this file and just have access to the
        :class:`~propka.molecular_container.MolecularContainer` object. If so,
        then pass ``False`` to ``write_pka``::

            mol = propka.run.single("protein.pdb", write_pka=False)

        In some cases, one may also want to pass a file-like (e.g.
        :class:`io.StringIO`) object instead of a file path as a string. In
        these cases the file-like object should be passed to the ``stream``
        argument and a string indicating the file type in the ``filename``
        argument; this string only has to look like a valid file name, it does
        not need to exist because the data are actually read from ``stream``.
        This approach is necessary because file-like objects do not usually
        have names, and propka uses the ``filename`` argument  to determine the
        input file type, and assigns the file name for the
        :class:`~propka.molecular_container.MolecularContainer` object::

            mol = propka.run.single('input.pdb', stream=string_io_file)

        In this case, a PDB file-like object was passed as `string_io_file`.
        The resultant pKa file will be written out as `input.pka`.

    Notes:
        * Only a single input structure file will be processed, defined by
          ``filename`` (and ``stream`` if passing a file-like object). Any
          additional files passed via the `-f` or `--file` flag to optargs will
          be ignored.


    .. seealso::

        :func:`propka.input.read_molecule_file`


    .. versionchanged:: 3.4.0
       Removed ability to write out PROPKA input files.
    """
    filename = str(filename)
    # Deal with input optarg options
    optargs = tuple(optargs)
    optargs += (filename,)
    options = loadOptions(optargs)

    parameters = read_parameter_file(options.parameters, Parameters())

    # Only filename present should be the one passed via the arguments
    # Anything else will probably have been passed using optargs' `-f` flag.
    ignored_list = [i for i in options.filenames if i != filename]
    if ignored_list:
        _LOGGER.warning(f"Ignoring extra filenames passed: {ignored_list}")
    options.filenames = [filename]

    my_molecule = MolecularContainer(parameters, options)
    my_molecule = read_molecule_file(filename, my_molecule, stream=stream)
    my_molecule.calculate_pka()

    # write outputs
    if write_pka:
        my_molecule.write_pka()

    return my_molecule
