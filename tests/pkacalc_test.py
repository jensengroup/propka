"""Reproduce tests from original PROPKA."""
import os
import re
from subprocess import call
import sys
import unittest
import logging
import propka.lib
import propka.molecular_container


# This error tolerance was chosen to make Ubuntu 18.04 pass under Windows 
# Subsystem Linux.
ACCEPTABLE_ERROR = 0.001


def run_propka(args):
    """Run PROPKA.

    Args:
        args:  command-line arguments for PROPKA
    """
    options = propka.lib.loadOptions()
    raise NotImplementedError("Need to incorporated command-line arguments")
    pdbfiles = options.filenames

    for pdbfile in pdbfiles:
        my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
        my_molecule.calculate_pka()
        my_molecule.write_pka()


# Setting this up as a direct translation of the original runtest.py script
# that will be run as part of 'python setup.py test'. This takes on the
# order of 10s to run.

class SystemTest(unittest.TestCase):
    """
    Run the program and compare against reference results.
    """
    def test_pka_calc(self):
        pdbs = ['1FTJ-Chain-A',
                '1HPX',
                '4DFR']

        test_dir = os.path.dirname(__file__)
        base_dir = os.path.dirname(test_dir)

        executable = os.path.join(base_dir, "scripts", "propka31.py")

        env = { "PYTHONPATH" : base_dir }

        for pdb in pdbs:
            input_filename = os.path.join(test_dir, "pdb", pdb + ".pdb")
            output_filename = os.path.join(test_dir, pdb + ".out")

            output_file = open(output_filename, "w")
            call([sys.executable, executable, input_filename],
                    stdout=output_file, env=env)
            output_file.close()

            # Check pka predictions.
            ref = open(os.path.join(test_dir, "results", pdb + ".dat"))
            output = open(output_filename)

            atpka = False
            errors = []
            for line in output:
                if not atpka:
                    # Start testing pka values.
                    if "model-pKa" in line:
                        atpka = True
                    continue

                m = re.search('([0-9]+\.[0-9]+)', line)
                if not m:
                    break

                expected_value = float(ref.readline())
                value = float(m.group(0))
                value_error = (value-expected_value)/expected_value

                if abs(value_error) > ACCEPTABLE_ERROR:
                    logging.error(value_error)
                    identity = line[:m.start()].strip()
                    errors.append("%12s  %8.2f    %8.2f" %
                            (identity, expected_value, value))

            os.remove("%s.pka" % pdb)
            os.remove("%s.propka_input" % pdb)

            ref.close()
            output.close()

            if errors:
                error_header = "       Group  Expected  Calculated\n"
                self.fail("Unexpected pKa values:\n" + error_header +
                          "\n".join(errors))
