# PROPKA 3.1

PROPKA predicts the pKa values of ionizable groups in proteins
(version 3.0) and protein-ligand complexes (version 3.1)
based on the 3D structure.

For proteins without ligands both version should produce the same result.

The method is described in the following papers, which you should cite
in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.

See [propka.ki.ku.dk](http://propka.ki.ku.dk/) for the PROPKA web server,
using the [tutorial](http://propka.ki.ku.dk/~luca/wiki/index.php/PROPKA_3.1_Tutorial).


## Installation

Clone repository or unpack the tar ball and install with
[setuptools](http://pythonhosted.org/setuptools/index.html) (note: if
you don't have setuptools installed you will need an internet
connection so that the installation procedure can download the
required files):

    cd propka-3.1
    python setup.py install --user

This will install the `propka31` script in your executable directory,
as configured for setuptools, for instance `~/.local/bin`. You can
change the bin directory with the `--install-scripts` option. For
example, in order to install in my `bin` directory in my home
directory:

    python setup.py install --user --install-scripts ~/bin


## Requirements

* Python 2.7 or higher or Python 3.1 or higher 

## Getting started

1. Clone the code from GitHub
2. `python setup.py install --user`
2. Run `propka31` with a .pdb file (see Examples)

## Examples

Calculate using pdb file

    propka31 1hpx.pdb


## Testing (for developers)

Please run `Tests/runtest.py/` after changes before pushing commits.

## References / Citations

Please cite these references in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.



