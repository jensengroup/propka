*********
Changelog
*********

3.5.1 (2023-12-31)
==================

Changes
-------

* Add more static typing and associated tests (`#177 <https://github.com/jensengroup/propka/pull/177>`_ and `#172 <https://github.com/jensengroup/propka/pull/172>`_)

* Refactor deprecated ``pkg_resources`` usage (`#176 <https://github.com/jensengroup/propka/pull/176>`_)

* Increase number of valence electrons for ligand atoms (`#169 <https://github.com/jensengroup/propka/pull/169>`_ and `#170 <https://github.com/jensengroup/propka/pull/170>_`)

Fixes
-----

* Fix topping up from all conformations, not just first (`#167 <https://github.com/jensengroup/propka/pull/167>`_)

3.5.0 (2023-02-13)
==================

Changes
-------

* Remove support for Python 3.6 and 3.7; add support for up to Python 3.11
  (`#154 <https://github.com/jensengroup/propka/issues/154>`_ and
  `#150 <https://github.com/jensengroup/propka/pull/150>`_)

* Add context manager for ``open_file_for_reading``
  (`#133 <https://github.com/jensengroup/propka/pull/133>`_)

Fixes
-----

* Fix precision of ``MolecularContainer.get_pi()``
  (`#148 <https://github.com/jensengroup/propka/pull/148>`_)

* Rename vanadium from ``Va`` to ``V``
  (`#141 <https://github.com/jensengroup/propka/pull/141>`_)

* Fix rounding issues in folding profile reporting
  (`#124 <https://github.com/jensengroup/propka/pull/124>`_)

* Fix malfunctioning ``-g -w`` command line options
  (`#124 <https://github.com/jensengroup/propka/pull/124>`_)


v3.4.0 (2020-12-19)
===================

Changes
-------

* Removed PROPKA input support and argument ``--generate-propka-input``
  (`#99 <https://github.com/jensengroup/propka/issues/99>`_)

* Add Python 3.9 support to continuous integration.
  (`#101 <https://github.com/jensengroup/propka/issues/101>`_)

* Removed logging abstraction from code to facilitate debugging and reduce code bloat.
  (`#108 <https://github.com/jensengroup/propka/issues/108>`_)


Fixes
-----

* Fixed bug that raised exception when missing amide nitrogen or oxygen.
  (`#17 <https://github.com/jensengroup/propka/issues/17>`_)

* ``propka --version`` now shows the program version and exits. Previously this option took a version argument to specify the sub-version of propka.
  However, this was non-functional at least since 2012.
  (`#89 <https://github.com/jensengroup/propka/issues/89>`_)

* Fix pI reporting in last line of :file:`.pka` file.
  (`<https://github.com/jensengroup/propka/pull/91>`_)

* Report correct version in :file:`.pka` file header.
  (`<https://github.com/jensengroup/propka/pull/92>`_)

* Fix handling of multi-model PDB without MODEL 1 entry.
  (`<https://github.com/jensengroup/propka/issues/96>`_)

* Fixed bug and sped up algorithm for identifying bonds via bounding boxes.
  (`#97 <https://github.com/jensengroup/propka/issues/97>`_, `#110 <https://github.com/jensengroup/propka/pull/110>`_)

* Fixed bug in ``propka --display-coupled-residues`` that crashed the program.
  (`#105 <https://github.com/jensengroup/propka/issues/105>`_)


v3.3.0 (2020-07-18)
===================

Additions
---------

* Add Sphinx documentation on `readthedocs.io <https://propka.readthedocs.io>`_
  (`#69 <https://github.com/jensengroup/propka/issues/69>`_, `#76 <https://github.com/jensengroup/propka/pull/76>`_, `#79 <https://github.com/jensengroup/propka/pull/79>`_)

Changes
-------

* Updated :func:`read_molecule_file` to accept file-like objects.
  (`#83 <https://github.com/jensengroup/propka/issues/83>`_)

* Use `versioneer <https://github.com/python-versioneer/python-versioneer>`_ for version management.
  (`#87 <https://github.com/jensengroup/propka/issues/87>`_)

* Add `code coverage <http://codecov.io>`_ to continuous integration pipeline.
  (`#62 <https://github.com/jensengroup/propka/pull/62>`_, `#71 <https://github.com/jensengroup/propka/pull/71>`_, `#76 <https://github.com/jensengroup/propka/pull/76>`_)

Fixes
-----

* Bundle required JSON files with package.
  (`#48 <https://github.com/jensengroup/propka/issues/48>`_)

* Fixed :class:`KeyError` bug in :func:`read_parameter_file`.
  (`#65 <https://github.com/jensengroup/propka/pull/65>`_)

* Update links to web server.
  (`#80 <https://github.com/jensengroup/propka/pull/80>`_)

* Fixed PDB reading for PROPKA "single" runs.
  (`#82 <https://github.com/jensengroup/propka/issues/82>`_)


v3.2.0 (2020-06-19)
===================

Additions
---------

* Significantly expanded testing framework.
  (`#30 <https://github.com/jensengroup/propka/pull/30>`_, `#36 <https://github.com/jensengroup/propka/pull/36>`_, `#37 <https://github.com/jensengroup/propka/pull/37>`_)

Changes
-------

* Improved ability to use PROPKA as a module in other Python scripts.
  (`#8 <https://github.com/jensengroup/propka/pull/8>`_)

* Improved output via :mod:`logging`.
  (`#11 <https://github.com/jensengroup/propka/pull/11>`_, `#12 <https://github.com/jensengroup/propka/pull/12>`_)

* Replaced data/parameter pickle file with human-readable JSON.
  (`#29 <https://github.com/jensengroup/propka/pull/29>`_)

* Significant delinting and formatting standardization against PEP8.
  (`#33 <https://github.com/jensengroup/propka/pull/33>`_, `#40 <https://github.com/jensengroup/propka/pull/40>`_)

* Improved package documentation.
  (`#41 <https://github.com/jensengroup/propka/pull/41>`_, `#61 <https://github.com/jensengroup/propka/pull/61>`_)

* Significant package refactoring.
  (`#46 <https://github.com/jensengroup/propka/issues/46>`_, `#47 <https://github.com/jensengroup/propka/pull/47>`_, `#59 <https://github.com/jensengroup/propka/pull/59>`_)

* Simplify module import structure.
  (`#49 <https://github.com/jensengroup/propka/issues/49>`_, `#61 <https://github.com/jensengroup/propka/pull/61>`_)

* Improved tempfile handling.
  (`#61 <https://github.com/jensengroup/propka/pull/61>`_)

v3.1.0
======

*Archaeologists wanted* to help us document the history of the code in versions 3.1.0 and earlier.