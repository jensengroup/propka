#!/usr/bin/env python

""" Run test for test pdbs """

from __future__ import division
from __future__ import print_function

from subprocess import call
import os, re
import sys

if __name__ == "__main__":
    # A list of input structures and command-line arguments to be passed in
    # to PROPKA for each:
    pdbs = [('1FTJ-Chain-A', []),
            ('1HPX', []),
            ('4DFR', []),
            ('3SGB', []),
            ('3SGB-subset', ['--titrate_only', 'E:17,E:18,E:19,E:29,E:44,E:45,E:46,E:118,E:119,E:120,E:139']),
            ('1HPX-info_warning', ['--log-level=30']),
           ]

    for pdb, args in pdbs:
        print('')
        print('RUNNING '+pdb)

        # Run pka calculation
        fh = open(pdb + '.out', 'w')
        cmd = [sys.executable, '../scripts/propka31.py','pdb/%s.pdb' % pdb] + args
        ret = call(cmd, stdout=fh, stderr=fh)
        if ret != 0:
            print(" ERR:")
            print(" Failed to execute PROPKA on %s" % pdb)
            print(" See: %s.out" % pdb)
            sys.exit(1)

        # Test pka predictions
        result_file = 'results/%s.dat' % pdb
        if not os.path.isfile(result_file):
            print(" ERR:")
            print(" file not found: %s" % result_file)
            sys.exit(1)

        result = open(result_file,'r')
        atpka = False
        for line in open(pdb+'.pka', 'r').readlines():
            if not atpka:
                if "model-pKa" in line:
                    # test pka
                    atpka = True
                    continue
                else:
                    continue
            if "-" in line:
                # done testing
                atpka = False
                continue

            expected_value = float(result.readline())
            m = re.search('([0-9]+\.[0-9]+)', line)
            value = float(m.group(0))

            if value != expected_value:
                print(" ERR:")
                print(line)
                print(" %s should be: %s" % (value, expected_value))
                sys.exit(1)

