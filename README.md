# PROPKA 3

PROPKA predicts the pKa values of ionizable groups in proteins
(version 3.0) and protein-ligand complexes (version 3.1 and later)
based on the 3D structure.

For proteins without ligands, both version should produce the same result.

The method is described in the following papers, which you should cite
in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295. doi:[10.1021/ct200133y](https://doi.org/10.1021/ct200133y)

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537. doi:[10.1021/ct100578z](https://doi.org/10.1021/ct100578z)

## PROPKA versions

The code in this repository is named _PROPKA 3_ and it is based on the original PROPKA 3.1 release (described in the papers above). It has undergone various changes, which is reflected in the version numbering. For instance, version 3.2 contains a number of API changes and code refactoring that introduce incompatibilities between the original 3.1 release and the more recent versions. In the future, we will increase the minor version number to indicate further changes to the code base (e.g., release 3.3 or 3.4). The major release number is not expected to change unless major changes to the underlying algorithms are implemented.

## Requirements

PROPKA 3 requires Python 3.6 or higher.  Additional requirements are specified in the `requirements.txt` file and automatically satisfied when installing with [pip](https://pip.pypa.io).

## Installation

### PIP-based installation

The easiest way to install PROPKA is via the [PyPI archive](https://pypi.org/project/PROPKA/) with the command

    pip install propka

This installation will install the `propka` Python module and the `propka3` executable script.
As always, a virtual environment (e.g., via [virtualenv](https://pypi.org/project/virtualenv/)) is recommended when installing packages.

### Source-based installation

The source code can be installed by cloning the repository or unpacking from a source code archive and running

    pip install .

in the source directory.
For the purposes of testing or development, you may prefer to install PROPKA as an editable module via PIP by running

    pip install -e .

in the source directory.

## Getting started

PROPKA can be used either as a module or via the installed script; i.e., either

    propka3

or

    python -m propka

works for invoking PROPKA.

A brief list of available options can be obtained by running PROPKA with no options:

    propka3

A longer list of options and descriptions is available using the `--help` option:

    propka3 --help

Most users run PROPKA by invoking the program with a PDB file as its argument; e.g.,

    propka3 1hpx.pdb

## Testing (for developers)

Please see [`tests/README.md`](tests/README.md) for testing instructions.
Please run these tests after making changes to the code and _before_ pushing commits.


## References / Citations

Please cite these references in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.
