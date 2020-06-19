.. -*- coding: utf-8 -*-

==============
 Installation
==============


PROPKA 3 requires Python 3.6 or higher. Additional requirements are
specified in the :file:`requirements.txt` file and automatically satisfied
when installing with pip_.


``pip``-based installation
==========================

The easiest way to install a release of PROPKA 3 is from the `PyPI archive`_ with the command

.. code-block:: bash

   pip install --upgrade propka

   
This installation will install the :mod:`propka` Python module and the
:command:`propka3` executable script. As always, a virtual environment (e.g., via
`virtualenv`_) is recommended when installing packages.

Source-based installation
=========================

The source code can be installed by cloning the `repository`_ or
unpacking from a source code archive and running

.. code-block:: bash

   pip install .

in the source directory.

For the purposes of testing or development, you may prefer to install
PROPKA as an editable module via pip_ by running

.. code-block:: bash

   pip install -e .

in the source directory.


.. _pip: https://pip.pypa.io/
.. _PyPI archive: https://pypi.org/project/PROPKA/
.. _virtualenv: https://pypi.org/project/virtualenv/
.. _repository: https://github.com/jensengroup/propka
