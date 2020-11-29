.. -*- coding: utf-8 -*-

.. _propka3-command:

============================
 :program:`propka3` command
============================

PROPKA predicts the pKa values of ionizable groups in proteins and
protein-ligand complexes based in the 3D structure. The
:program:`propka3` command has the following options::
	   
  propka3 [-h] [-f FILENAMES] [-r REFERENCE] [-c CHAINS] [-i TITRATE_ONLY] [-t THERMOPHILES] [-a ALIGNMENT] [-m MUTATIONS]
         [-v VERSION_LABEL] [-p PARAMETERS] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [-o PH] [-w WINDOW WINDOW WINDOW]
         [-g GRID GRID GRID] [--mutator MUTATOR] [--mutator-option MUTATOR_OPTIONS] [-d] [-l] [-k] [-q] [--protonate-all]
         input_pdb


.. program:: propka3
	     
.. option:: input_pdb

            read data from file <input_pdb>

.. option::   -h, --help

	      show this help message and exit
	      
.. option::  -f FILENAMES, --file FILENAMES
	     
             read data from <filename>, i.e. <filename> is added to
	     arguments (default: [])
	     
.. option::  -r REFERENCE, --reference REFERENCE
	     
             setting which reference to use for stability calculations
	     [neutral/low-pH] (default: neutral)
	     
.. option::  -c CHAINS, --chain CHAINS
	     
             creating the protein with only a specified chain. Specify
	     " " for chains without ID [all] (default: None)
	     
.. option::  -i TITRATE_ONLY, --titrate_only TITRATE_ONLY
	     
             Treat only the specified residues as titratable. Value
	     should be a comma-separated list of "chain:resnum"
	     values; for  example: ``-i "A:10,A:11"`` (default: None)
	     
.. option::  -t THERMOPHILES, --thermophile THERMOPHILES
	     
             defining a thermophile filename; usually used in
	     'alignment-mutations' (default: None)
	     
.. option::  -a ALIGNMENT, --alignment ALIGNMENT
	     
             alignment file connecting <filename> and <thermophile>
	     [<thermophile>.pir] (default: None)
	     
.. option::  -m MUTATIONS, --mutation MUTATIONS
	     
              specifying mutation labels which is used to modify
	      <filename> according to, e.g. N25R/N181D (default: None)
			
.. option::  -v VERSION_LABEL, --version VERSION_LABEL
	     
             specifying the sub-version of propka [Jan15/Dec19]
	     (default: Jan15)
	     
.. option::  -p PARAMETERS, --parameters PARAMETERS
	     
             set the parameter file (default:
	     <installation_directory>/propka/propka/propka.cfg)
	     
.. option::  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
	     
             logging level verbosity (default: INFO)
	     
.. option::  -o PH, --pH PH

	     setting pH-value used in e.g. stability calculations
	     (default: 7.0)
	     
.. option::  -w WINDOW WINDOW WINDOW, --window WINDOW WINDOW WINDOW
	     
             setting the pH-window to show e.g. stability profiles
	     (default: (0.0, 14.0, 1.0)) 

.. option:: -g GRID GRID GRID, --grid GRID GRID GRID
	    
            setting the pH-grid to calculate e.g. stability related
	    properties (default: (0.0, 14.0, 0.1))
	    
.. option::  --mutator MUTATOR

	     setting approach for mutating <filename>
	     [alignment/scwrl/jackal] (default: None)
	     
.. option::  --mutator-option MUTATOR_OPTIONS
	     
             setting property for mutator [e.g. type="side-chain"]
	     (default: None)
	     
.. option::  -d, --display-coupled-residues
	     
             Displays alternative pKa values due to coupling of
	     titratable groups (default: False)
	     
.. option::  -l, --reuse-ligand-mol2-files
	     
             Reuses the ligand mol2 files allowing the user to alter
	     ligand bond orders (default: False)
	     
.. option::  -k, --keep-protons

	     Keep protons in input file (default: False)
	     
.. option::  -q, --quiet
	     
             suppress non-warning messages (default: None)
	     
.. option::  --protonate-all

	     Protonate all atoms (will not influence pKa calculation)
	     (default: False)
	     
