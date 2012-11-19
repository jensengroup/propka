# PROPKA 3.1

PROPKA predicts the pKa values of ionizable groups in
proteins (version 3.0) and
protein-ligand complexes (version 3.1)
based in the 3D structure.

For proteins without ligands both version should produce the same result.

The method is described in the following papers, which you should cite
in publications:

* Søndergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Søndergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.

See [propka.ki.ku.dk](http://propka.ki.ku.dk/) for the PROPKA web server,
using the [tutorial](http://propka.ki.ku.dk/~luca/wiki/index.php/PROPKA_3.1_Tutorial).

## Installation

No installation needed. Just clone and run.

## Requirements

* Python 3.1 or higher

## Getting started

1. Clone the code from GitHub
2. Run 'propka.py' with a .pdb file (see Examples)

## Examples

Calculate using pdb file

    ./propka.py 1hpx.pdb

If for some reason your setup with python3.1+ is
not located in '/usr/bin/python3', run the script

    python3.2 propka.py 1hpx.pdb

## Testing (for developers)

Please run `Tests/runtest.py/` after changes before pushing commits.

## References / Citations

Please cite these references in publications:

* Søndergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Søndergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.



