# PROPKA 3.2

PROPKA predicts the pKa values of ionizable groups in proteins
(version 3.0) and protein-ligand complexes (version 3.2 and later)
based on the 3D structure.

For proteins without ligands, both version should produce the same result.

The method is described in the following papers, which you should cite
in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295. doi:[10.1021/ct200133y](https://doi.org/10.1021/ct200133y)

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537. doi:[10.1021/ct100578z](https://doi.org/10.1021/ct100578z)

## Requirements

PROPKA 3 requires Python 3.5 or higher.  Additional requirements are specified in the `requirements.txt` file and automatically satisfied when installing with [pip](https://pip.pypa.io).

## Installation

### PIP-based installation

The easiest way to install PROPKA is via the [PyPI archive](https://pypi.org/project/PROPKA/) with the command

    pip install propka

As always, a virtual environment (e.g., via [virtualenv](https://pypi.org/project/virtualenv/)) is recommended when installing packages.

For the purposes of testing or development, you may prefer to install PROPKA as
an editable module via PIP by running

    pip install -e .

### Source-based installation

The source code can be installed by cloning the repository or unpacking from a source code archive and running

    pip install .

in the source directory.

Installation is also possible with 
[setuptools](http://pythonhosted.org/setuptools/index.html) which offers additional customization options; e.g.

    python setup.py install --user

will install the `propka3` script in your executable directory.

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
