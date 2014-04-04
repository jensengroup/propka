# PROPKA 3.1
#
#
# setuptools installation of  PROPKA 3.1

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages

VERSION = "3.1"

setup(name="PROPKA",
      version=VERSION,
      description="Heuristic pKa calculations with ligands",
      long_description="""
PROPKA predicts the pKa values of ionizable groups in proteins (version 3.0) and
protein-ligand complexes (version 3.1) based on the 3D structure.

For proteins without ligands both version should produce the same result.

The method is described in the following papers, which you should cite
in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan
  H. Jensen. "Improved Treatment of Ligands and Coupling Effects in
  Empirical Calculation and Rationalization of pKa Values." Journal of
  Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan
  H. Jensen. "PROPKA3: consistent treatment of internal and surface
  residues in empirical pKa predictions." Journal of Chemical Theory
  and Computation 7, no. 2 (2011): 525-537.

See http://propka.ki.ku.dk/ for the PROPKA web server,
using the tutorial http://propka.ki.ku.dk/~luca/wiki/index.php/PROPKA_3.1_Tutorial .
""",
      author="Jan H. Jensen",
      author_email="jhjensen@chem.ku.dk",
      license="",
      url="http://propka.ki.ku.dk/",
      keywords="science",
      packages=find_packages(exclude=['scripts']),
      package_data = {'propka': ['*.dat', '*.cfg']},
      #scripts = ["scripts/propka31.py"],  # use entry point below
      entry_points = {
        'console_scripts': [
            'propka31 = propka.run:main',
            ],
        },
      zip_safe=True,
      test_suite="Tests",
)
