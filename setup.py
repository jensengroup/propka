#! /usr/bin/python
# PROPKA 3

from setuptools import setup, find_packages

VERSION = "3.2.0"

setup(
    name="propka",
    version=VERSION,
    description="Heuristic pKa calculations with ligands", long_description="""
PROPKA predicts the pKa values of ionizable groups in proteins (version 3.0) and
protein-ligand complexes (version 3.1 and later) based on the 3D structure.

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

See http://propka.org/ for the PROPKA web server.
""",
    author="Jan H. Jensen",
    author_email="jhjensen@chem.ku.dk",
    maintainer="Nathan Baker",
    maintainer_email="nathanandrewbaker@gmail.com",
    license="LGPL v2.1",
    url="http://propka.org",
    keywords="science",
    classifiers=[
        'Development Status :: 6 - Mature',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Operating System :: POSIX', 
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        ],
    packages=find_packages(exclude=['scripts']),
    package_data={'propka': ['*.dat', '*.cfg', '*.json']},
    entry_points={'console_scripts': ['propka3 = propka.run:main', ]},
    zip_safe=True,
    python_requires='>=3.6',
    tests_require=["pandas", "numpy", "pytest"],
    test_suite="tests",
    )
