.. -*- coding: utf-8 -*-



.. |pKa| replace:: :math:`\text{p}K_\text{a}`

==================
 Quickstart Guide
==================

PROPKA can be used either via the installed script :program:`propka3`
or as a Python module. When using the :ref:`propka3-command`, use

.. code-block:: bash

   propka3  FILENAME

As a module (:mod:`propka`), also provide the input filename

.. code-block:: bash
		
   python -m propka FILENAME

In both cases, additional options may be added, as described in more
detail for the :ref:`propka3-command`.


Predicting protein residue |pKa| values
=======================================

Most users run PROPKA by invoking the :program:`propka3` program with
a PDB file as its argument; e.g., for PDB 1HPX_ (HIV-1 protease
complexed with the inhibitor KNI-272)

.. code-block:: bash

   propka3 1hpx.pdb

In this example, |pKa| values of titratable protein residues and
titratable groups of the inhibitor KNI are calculated.
   
The output looks similar to the following (many lines omitted as
"..."). It is also contained in the output file :file:`1hpx.pka` that
is automatically written::


   propka3.2                                                                                    2020-06-19
   ...
   ...
   Found NAR group:  1530-  N1   900-KNI (B) [   7.907    1.459    5.427] N
   Found O3  group:  1531-  O1   900-KNI (B) [   5.235    3.791    9.082] O
   Found O2  group:  1532-  O3   900-KNI (B) [   3.327    4.297   11.852] O
   Found NAM group:  1533-  N2   900-KNI (B) [   3.955    2.384   10.893] N
   Found O2  group:  1539-  O6   900-KNI (B) [   3.758   -0.629   12.111] O
   Found NAM group:  1541-  N3   900-KNI (B) [   4.496    0.982   13.492] N
   Found O2  group:  1542-  O4   900-KNI (B) [   6.324   -1.234   17.045] O
   Found OH  group:  1548-  O2   900-KNI (B) [   4.949    0.934   16.427] O
   Found O2  group:  1559-  O5   900-KNI (B) [   6.746   -3.574   14.588] O
   Found NAM group:  1560-  N5   900-KNI (B) [   7.637   -4.575   16.403] N
   -------------------------------------------------------------------------------------------------------


   Calculating pKas for Conformation container 1A with 1878 atoms and 480 groups
   -------------------------------------------------------------------------------------------------------

   ---------  -----   ------   ---------------------    --------------    --------------    --------------
			       DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC    
    RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION  
   ---------  -----   ------   ---------   ---------    --------------    --------------    --------------

   ASP  25 A   5.07*  100 %    4.30  617   0.19    0   -0.85 KNI  O4 B   -0.63 GLY  27 A    0.07 ASP  29 A
   ASP  25 A                                           -0.85 KNI  O2 B   -0.09 ALA  28 A    0.00 XXX   0 X
   ASP  25 A                                           -0.84 ASP  25 B   -0.04 GLY  27 B    0.00 XXX   0 X

   ASP  29 A   3.11    50 %    1.20  420   0.13    0   -0.68 ARG  87 A    0.00 XXX   0 X   -0.04 LYS  45 A
   ASP  29 A                                           -0.28 ARG   8 B    0.00 XXX   0 X   -0.47 ARG  87 A
   ASP  29 A                                            0.00 XXX   0 X    0.00 XXX   0 X   -0.54 ARG   8 B

   ASP  30 A   4.62    59 %    1.30  446   0.00    0   -0.11 LYS  45 A    0.00 XXX   0 X   -0.07 ARG  87 A
   ASP  30 A                                            0.00 XXX   0 X    0.00 XXX   0 X   -0.01 ARG   8 B
   ASP  30 A                                            0.00 XXX   0 X    0.00 XXX   0 X    0.29 ASP  29 A
   ASP  30 A                                            0.00 XXX   0 X    0.00 XXX   0 X   -0.57 LYS  45 A

   ASP  60 A   2.55     0 %    0.41  249   0.00    0   -0.40 THR  74 A    0.00 XXX   0 X   -0.02 LYS  45 A
   ASP  60 A                                           -0.85 LYS  43 A    0.00 XXX   0 X   -0.38 LYS  43 A
   ...
   ...
   ...
   ARG  87 B  12.28    45 %   -1.40  407   0.00    0    0.77 ASP  29 B    0.00 XXX   0 X    0.10 ASP  30 B
   ARG  87 B                                            0.00 XXX   0 X    0.00 XXX   0 X   -0.19 ARG   8 A
   ARG  87 B                                            0.00 XXX   0 X    0.00 XXX   0 X    0.50 ASP  29 B

   N+    1 B   8.96     0 %   -0.39  235   0.00    0    0.85 C-   99 A    0.00 XXX   0 X    0.07 CYS  67 B
   N+    1 B                                            0.00 XXX   0 X    0.00 XXX   0 X    0.04 CYS  95 B
   N+    1 B                                            0.00 XXX   0 X    0.00 XXX   0 X    0.38 C-   99 A

   KNI  N1 B   4.60     0 %   -0.36  273   0.00    0    0.00 XXX   0 X    0.00 XXX   0 X   -0.03 ARG   8 A

   Coupled residues (marked *) were detected.Please rerun PropKa with the --display-coupled-residues 
   or -d option for detailed information.

   --------------------------------------------------------------------------------------------------------
   SUMMARY OF THIS PREDICTION
	  Group      pKa  model-pKa   ligand atom-type
      ASP  25 A     5.07       3.80                      
      ASP  29 A     3.11       3.80                      
      ASP  30 A     4.62       3.80                      
      ASP  60 A     2.55       3.80                      
      ASP  25 B     9.28       3.80                      
      ASP  29 B     1.78       3.80                      
      ASP  30 B     4.91       3.80                      
      ASP  60 B     2.13       3.80                      
      GLU  21 A     4.78       4.50                      
      GLU  34 A     3.93       4.50                      
      GLU  35 A     3.65       4.50                      
      GLU  65 A     3.89       4.50                      
      GLU  21 B     4.73       4.50                      
      GLU  34 B     3.36       4.50                      
      GLU  35 B     4.07       4.50                      
      GLU  65 B     3.70       4.50                      
      C-   99 A     2.08       3.20                      
      C-   99 B     2.11       3.20                      
      HIS  69 A     6.98       6.50                      
      HIS  69 B     7.11       6.50                      
      CYS  67 A     9.41       9.00                      
      CYS  95 A    11.68       9.00                      
      CYS  67 B     9.82       9.00                      
      CYS  95 B    11.61       9.00                      
      TYR  59 A     9.67      10.00                      
      TYR  59 B     9.54      10.00                      
      LYS  14 A    10.43      10.50                      
      LYS  20 A    10.32      10.50                      
      LYS  43 A    11.41      10.50                      
      LYS  45 A    10.54      10.50                      
      LYS  55 A    10.42      10.50                      
      LYS  70 A    10.92      10.50                      
      LYS  14 B    10.55      10.50                      
      LYS  20 B    11.01      10.50                      
      LYS  43 B    11.43      10.50                      
      LYS  45 B    10.47      10.50                      
      LYS  55 B    10.41      10.50                      
      LYS  70 B    11.07      10.50                      
      ARG   8 A    13.96      12.50                      
      ARG  41 A    12.41      12.50                      
      ARG  57 A    14.40      12.50                      
      ARG  87 A    12.35      12.50                      
      ARG   8 B    12.76      12.50                      
      ARG  41 B    12.42      12.50                      
      ARG  57 B    13.73      12.50                      
      ARG  87 B    12.28      12.50                      
      N+    1 A     8.96       8.00                      
      N+    1 B     8.96       8.00                      
      KNI  N1 B     4.60       5.00                NAR   

   Writing 1hpx.pka



Some of the important contents:

- The section *Calculating pKas for Conformation container 1A with 1878 atoms and
  480 groups* lists details on the calculations for all ionizable
  residues. It indicates the considerations that went into a |pKa|
  estimate such as hydrogen bonds and Coulomb interactions. It also
  indicates if there is potentially coupling between residues.
- Values with "XXX" placeholders are not calculated (but appear to
  maintain the formatting).
- The section *SUMMARY OF THIS PREDICTION* lists the predicted |pKa|
  for each residue together with the model |pKa| (the "default"
  value).
- Ligand values are labeled with the residue name of the ligand, in
  this case "KNI".  


.. todo::

   More detailed discussion of output.
  


Getting help
============

A brief list of available options can be obtained by running PROPKA
with no options. A longer list of options and descriptions is
available using the :option:`propka3 --help` option:

.. code-block:: bash

   propka3 --help

   
.. links
.. _1HPX: https://www.rcsb.org/structure/1HPX   
