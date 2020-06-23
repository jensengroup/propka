.. -*- coding: utf-8 -*-

===============
 API Reference
===============

The :program:`propka3` command provides a command-line interface to
PROPKA 3's functionality. It is built on classes and functions in the
:mod:`propka` module. The API of :mod:`propka` is documented here for
developers who might want to directly use the PROPKA 3 code.

.. Note::

   The API is still changing and there is currently no guarantee that
   it will remain stable between minor releases.



.. currentmodule:: propka

.. module:: propka		   


Data structures
===============

.. autosummary::
   :toctree: api
      
   atom
   bonds
   group   
   conformation_container
   molecular_container

   
I/O
===

.. autosummary::
   :toctree: api

   input
   lib
   output
   parameters
   hybrid36
   ligand_pka_values
   run
   version
   
   
Structure processing
====================

.. autosummary::
   :toctree: api

   protonate
   hydrogens	     
   ligand

	     
Calculations
============

.. autosummary::
   :toctree: api

   calculations
   coupled_groups
   determinant
   determinants
   energy
   iterative
   vector_algebra


