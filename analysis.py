#!/usr/bin/env python
"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of BondMatcher.

BondMatcher is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Some analysis of the network.
"""

import argparse
from collections import defaultdict

from src import files_io

parser = argparse.ArgumentParser(description='Simple chain analysis.')

parser.add_argument('-f', '--pdb', help='The path to the pdb file to read.', required=True)
parser.add_argument('-t', '--top', help='The path to the original topology file', required=True)
parser.add_argument('-c', '--config', help='The .par config file', required=True)

args = parser.parse_args()

if args.pdb.endswith('pdb'):
  pdb_file = files_io.PDBFile(args.pdb)
elif args.pdb.endswith('gro'):
  pdb_file = files_io.GROFile(args.pdb)
else:
  raise ValueError('Unrecognized file type, supported: .pdb, .gro')

print 'PDB file: %s\nTOP file: %s\nConfig file: %s' % (
    args.pdb, args.top, args.config
    )

# Reads the config.
execfile(args.config)

# Reads the input files
top_file = files_io.TopologyFile(args.top)

pdb_file.open()
top_file.open()

pdb_file.read()
top_file.read()


def calculate_linkage():
  """Calculate the linkage"""
  chain_name = raw_input('chain name: ')

  chain_pairs = set()
  for b in CHAIN_BONDS:  # pylint: disable=E0602
    chain_pairs.add((b[0][1], b[1][1]))
  chain_idx_bonds = defaultdict(list)

  number_of_loops = 0

  # Lets check the bonds
  for bond in top_file.bonds:
    ai = pdb_file.data[bond[0]]
    aj = pdb_file.data[bond[1]]
    chain_pair_1 = (ai.chain_name, aj.chain_name)
    chain_pair_2 = (aj.chain_name, ai.chain_name)
    if chain_pair_1 in chain_pairs or chain_pair_2 in chain_pairs:
      chain_idx, chain_aj_idx = (
          (ai.chain_idx, aj.chain_idx) if ai.chain_name == chain_name
          else (aj.chain_idx, ai.chain_idx)
          )
      if chain_aj_idx in chain_idx_bonds[chain_idx]:
        number_of_loops += 1
        print 'chain %s(%s) is already linked with %s(%s)' % (
            chain_idx, ai.atom_id, chain_aj_idx, aj.atom_id)
      chain_idx_bonds[chain_idx].append(chain_aj_idx)

  linkedge_count = defaultdict(int)
  for val in chain_idx_bonds.values():
    linkedge_count[len(val)] += 1

  if linkedge_count:
    print ''
    print 'number_of_links;number_of_%s_chains' % chain_name
    for count, value in linkedge_count.iteritems():
      print '%d;%d' % (count, value)
    print ''
  if linkedge_count:
    print 'Number of chains that formed loops: %d' % number_of_loops
    yn = raw_input('Do you want to save in file? [y/n]: ')
    if yn.lower().startswith('y'):
      f_path = raw_input('The output file: ')
      f_st = open(f_path, 'w')
      f_st.writelines('number_of_links;number_of_%s_chains\n' % chain_name)
      for count, value in linkedge_count.iteritems():
        f_st.writelines('%d;%d\n' % (count, value))
      f_st.close()


def number_of_bonds():
  """Calculates the number of bonds."""

  chain_pairs = set()
  for b in CHAIN_BONDS:  # pylint: disable=E0602
    chain_pairs.add((b[0][1], b[1][1]))

  nr_of_bonds = 0

  for bond in top_file.bonds:
    ai = pdb_file.data[bond[0]]
    aj = pdb_file.data[bond[1]]
    chain_pair_1 = (ai.chain_name, aj.chain_name)
    chain_pair_2 = (aj.chain_name, ai.chain_name)
    if chain_pair_1 in chain_pairs or chain_pair_2 in chain_pairs:
      nr_of_bonds += 1
  print 'Number of bonds between chain pairs: %d' % nr_of_bonds


def side_chains():
  pass


options = {
    1: ('Calculate the number of linked chains', calculate_linkage),
    2: ('Number of bonds', number_of_bonds),
    3: ('Side chains', side_chains),
    }

user_opt_id = None
while user_opt_id != 0:
  print 'Available options'
  print '=' * 17
  for opt_id, opt_desc in options.iteritems():
    print '%d - %s' % (opt_id, opt_desc[0])
  print '0 - Exit\n'
  user_input = raw_input('> ')
  user_opt_id = int(user_input)
  if user_opt_id != 0:
    options[int(user_opt_id)][1]()
