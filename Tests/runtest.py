#!/usr/bin/python3

""" Run test for test pdbs """

from __future__ import division
from __future__ import print_function

from subprocess import call
import os, re
import sys

pdbs = ['1FTJ-Chain-A',
        '1HPX',
        '4DFR']

for pdb in pdbs:
  print('')
  print('RUNNING '+pdb)

  # Run pka calculation
  call([sys.executable, '../propka.py','pdb/'+pdb+'.pdb'], stdout = open(pdb+'.out', 'w+'))

  # Test pka predictiona
  result = open('results/'+pdb+'.dat','r')
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

    r = float(result.readline())
    m = re.search('([0-9]+\.[0-9]+)', line)

    if(float(m.group(0)) != r):
      print(" ERR:")
      print(line)
      print(" "+"should be: "+str(r))






