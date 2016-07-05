#!/usr/bin/env python
"""
Copyright (C) 2014-2016 Jakub Krajniak (jkrajniak@gmail.com)

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

BondMatcher - main file
"""


import argparse
import sys
import time

from src import analysis
from src import files_io
from src import structures
from src import tools


# Define the arguments.
parser = argparse.ArgumentParser(description='Find bonds to join together.')
parser.add_argument('-f', '--pdb', help='The path to the pdb/gro file to read.', required=True)
parser.add_argument('-t', '--top', help='The path to the original topology file', required=True)
parser.add_argument('-o', '--output', help='The path to the output topology file', required=True)
parser.add_argument('-c', '--config', help='The name of the config class', required=True)
parser.add_argument('-v', '--verbose', help='Verbose mode', default=False, action='store_true')
parser.add_argument('--search_cutoff', help='Search for optimal cutoff', default=False,
                    action='store_true')
parser.add_argument('-d', '--density', help='Define the density [range: 0.0-1.0]', required=False,
                    type=float)
parser.add_argument('-co', '--cutoff', help='Init cutoff [A]', required=False, type=float)
parser.add_argument('-cd', '--cutoff_dt', help='Cutoff change step [A]', type=float)
parser.add_argument('-st', '--stats_file', help='Stats file', required=False)
parser.add_argument('-lf', '--log_file', help='The log file, by default run_log.log',
                    default='run_log.log')

args = parser.parse_args()

log_file = files_io.prepare_path(args.log_file)

logger = tools.Logger('', log_file)
logger.enable()

# Save the execution parameters in the log file.
logger.write_to_file('\nProgram executed with the parameters: %s\n' % ' '.join(sys.argv))

# Reads the config.
local_settings = structures.BondSettings(args.config)

if args.search_cutoff:
  if args.density or local_settings.EXPECTED_DENSITY:
    if not args.cutoff and not local_settings.INITIAL_CUTOFF:
      print 'Please define the init cut-off'
      sys.exit(1)

    if not args.cutoff_dt and not local_settings.CUTOFF_STEP:
      print 'Please define the cut-off step'
      sys.exit(1)

# Times
time0 = time.time()

if args.pdb.endswith('pdb'):
  pdb_file = files_io.PDBFile(args.pdb)
elif args.pdb.endswith('gro'):
  pdb_file = files_io.GROFile(args.pdb)
else:
  raise ValueError('Unrecognized file type, supported: .pdb, .gro')

print 'PDB/GRO file: %s\nTOP file: %s\nConfig file: %s\nOutput file: %s' % (
    args.pdb, args.top, args.config, args.output
    )

# Presave in the log file the config file, for the purpose of documenting the
# science.
logger.write_to_file('\n\n===  CONFIG FILE: [%s]  ===\n' % args.config)
for conf_l in open(args.config).readlines():
  if conf_l and not conf_l.startswith('#'):
    logger.write_to_file(conf_l)

logger.write_to_file('\n\n===   END CONFIG FILE  ===\n')

top_file = files_io.TopologyFile(args.top)

pdb_file.open()
top_file.open()

pdb_file.read()
top_file.read()

local_settings.BOX = pdb_file.box
print 'Simulation box: %s' % pdb_file.box

time1 = time.time()

# Preprocess atoms
valid_atoms, number_of_bonds, valid_bonds = tools.prepare_atoms(
    pdb_file.data,
    local_settings,
    top_file
    )


chain_ids_in_the_bonds, chain_types_in_bonds = tools.find_chains(
    top_file,
    pdb_file,
    local_settings
    )

print '\nBasic properties: '
print '- number of all atoms: %d\n- number of active sites: %d\n- number of formed bonds: %d' % (
    len(pdb_file.data), len(valid_atoms), number_of_bonds
    )
density = number_of_bonds / float(local_settings.TOTAL_NUMBER_OF_CROSSLINKS)
print '- density: %.4f' % density
print '- number of chains in network: %d' % len(chain_ids_in_the_bonds)
for chain_name, val in chain_types_in_bonds.iteritems():
  print '  - %s: %d' % (chain_name, len(val))

# Get the cutoff
if args.cutoff:
  cutoff = args.cutoff
elif args.search_cutoff:
  cutoff = local_settings.INITIAL_CUTOFF
else:
  cutoff = tools.get_cutoff(density, local_settings.CUTOFF_DISTANCE, 4.0)

if args.density and not args.search_cutoff:
  raise ValueError('Defining expected density without searching for cutoff does not make sense.')

if args.density:
  expected_density = args.density
  if expected_density < 0.0 or expected_density > 1.0:
    raise ValueError('Density %.2f is outside the valid range [0.0;1.0]' % expected_density)
elif args.search_cutoff:
  expected_density = local_settings.EXPECTED_DENSITY
else:
  expected_density = None

if expected_density:
  if density >= expected_density:
    print 'The current density %.2f is higher/equal that the defined %.2f' % (
        density, expected_density)
    sys.exit(1)

if args.cutoff_dt:
  cutoff_dt = args.cutoff_dt
elif args.search_cutoff:
  cutoff_dt = local_settings.CUTOFF_STEP
else:
  cutoff_dt = None

time2 = time.time()

if args.search_cutoff:
  print 'Looking for density %.2f %% with initial cutoff %.2f and cutoff step %.4f' % (
      expected_density * 100.0, cutoff, cutoff_dt
      )

# Presave in the log file, the execution config that can depend on the arguments
# from the command line and those save in the file.
logger.write_to_file(
"""
EXECUTION PARAMETERS:
- performing search cutoff: %s
- expected crosslink density: %s
- initial cutoff distance: %s
- cutoff step: %s
""" % (args.search_cutoff, expected_density, cutoff, cutoff_dt))

print '\nStart searching\n'
# Track the number of new bonds in respect to the cutoff
cutoff_stats = {}
search_for_bonds = True
step_ids = 0

# Analysis tool
analysis_tool = analysis.Analysis(
    None,
    local_settings,
    pdb_file.chains
    )

analysis_tool.add_bonds(valid_bonds)
logger.disable()

while search_for_bonds:
  if cutoff not in cutoff_stats:
    cutoff_stats[cutoff] = {}

  time_1 = time.time()
  print '%d: Looking for the bonds with cutoff with the atoms (%d): %.2f' % (
      step_ids, len(valid_atoms), cutoff
      )
  new_bonds = analysis_tool.find_atom_pairs(valid_atoms, cutoff)

  # Update the set of valid atoms, based on the the degree.
  valid_atoms = [x for x in valid_atoms if x.degree < x.max_degree]

  number_of_new_bonds = 0
  if args.verbose:
    print 'New bonds: '
  if new_bonds:
    for atom_pair, bonds in new_bonds.iteritems():
      bond_config = local_settings.BOND_CONFIG[atom_pair]
      for atom_i, atom_j in bonds:
        if (atom_i.atom_id, atom_j.atom_id) not in top_file.bonds:
          top_file.new_data['bonds'][(atom_i.atom_id, atom_j.atom_id)] = bond_config
          number_of_new_bonds += 1
        atom_i.neighbours.append(atom_j)
        atom_j.neighbours.append(atom_i)
        chain_ids_in_the_bonds.add(atom_j.chain_idx)
        chain_ids_in_the_bonds.add(atom_i.chain_idx)
        chain_types_in_bonds[atom_i.chain_name].add(atom_i.chain_idx)
        chain_types_in_bonds[atom_j.chain_name].add(atom_j.chain_idx)

        if args.verbose:
          print '%2s(id: %5s, %s, %s) - %2s(id: %5s, %s, %s)' % (
              atom_i.name, atom_i.atom_id, atom_i.chain_name, atom_i.position,
              atom_j.name, atom_j.atom_id, atom_j.chain_name, atom_j.position)
  # Prepare new entries for angles and dihedrals. Only if there are some new bonds.
  if number_of_new_bonds > 0:
    analysis_tool.bonds = new_bonds
    analysis_tool.add_bonds(new_bonds)

    new_angles = analysis_tool.find_angles()
    print 'Found %d new angles' % len(new_angles)
    if args.verbose:
      for ang_def in new_angles:
        a1, a2, a3 = [top_file.atoms[x] for x in ang_def[:3]]
        print '(%s[%s], %s[%s]) - (%s[%s], %s[%s]) - (%s[%s], %s[%s])' % (
          a1.atom_id, a1.chain_idx, a1.name, a1.chain_name,
          a2.atom_id, a2.chain_idx, a2.name, a2.chain_name,
          a3.atom_id, a3.chain_idx, a3.name, a3.chain_name
          )
    top_file.new_data['angles'].extend(new_angles)
    new_dihedrals = analysis_tool.find_dihedrals()
    print 'Found %d new dihedrals angles' % len(new_dihedrals)
    top_file.new_data['dihedrals'].extend(new_dihedrals)
    if args.verbose:
      for dih_def in new_dihedrals:
        a1, a2, a3, a4 = [top_file.atoms[x] for x in dih_def[:4]]
        print '%s[%s] %s[%s] %s[%s] %s[%s]' % (
          a1.name, a1.chain_name,
          a2.name, a2.chain_name,
          a3.name, a3.chain_name,
          a4.name, a4.chain_name
          )


  # Save some stats
  cutoff_stats[cutoff]['number_of_bonds'] = number_of_bonds
  number_of_bonds += number_of_new_bonds
  density = float(number_of_bonds) / local_settings.TOTAL_NUMBER_OF_CROSSLINKS

  print ('%d: Total number of bonds: %d, '
         'number of new bonds: %d, '
         'current cross-link density %.2f% %, '
         'number of chains %s, '
         'total_number: %d'
        ) % (step_ids,
             number_of_bonds,
             number_of_new_bonds,
             density * 100.0,
             ', '.join(['%s:%d' % (k, len(v)) for (k, v) in chain_types_in_bonds.iteritems()]),
             len(chain_ids_in_the_bonds)
             )

  # Save some stats
  cutoff_stats[cutoff]['number_of_new_bonds'] = number_of_new_bonds
  cutoff_stats[cutoff]['density'] = density
  cutoff_stats[cutoff]['number_of_chains'] = len(chain_ids_in_the_bonds)
  cutoff_stats[cutoff].update({
      k: len(v) for (k, v) in chain_types_in_bonds.iteritems()
      })

  if args.search_cutoff:
    if density < expected_density:
      cutoff += cutoff_dt
    else:
      search_for_bonds = False
  else:
    search_for_bonds = False
  step_ids += 1

# Display some timing information
if args.verbose:
  print '\n= Time performance ='
  print 'Filtering bonds: %.4fs' % (time2 - time1)
  print 'Searching bonds: %.4fs' % (time.time() - time2)
  print 'Total execution time: %.4fs' % (time.time() - time0)

print '\nFinished!'
print 'Total number of bonds: %d' % number_of_bonds
print 'Cross-link density: %.2f %%' % (density * 100.0)
print 'Cutoff distance: %.4f' % cutoff


#
# Postprocess the topology
#
# Looking for the pairs
new_pairs = analysis.find_1_4_pairs(
    analysis_tool.all_bonds,
    local_settings
    )

# Filter new pairs
filtered_new_pairs = {
    key: value for key, value in new_pairs.iteritems()
    if (key[0], key[1]) not in top_file.pairs and (key[1], key[0]) not in top_file.pairs
    }

print 'Found in total %d of new pairs' % len(filtered_new_pairs)
top_file.new_data['pairs'] = filtered_new_pairs

#
# Write the topology file
#
top_file.write(args.output)

# Save some stats
if args.stats_file:
  header = [
      'number_of_bonds',
      'number_of_new_bonds',
      'density',
      'number_of_chains'
      ] + chain_types_in_bonds.keys()
  stats_file = open(args.stats_file, 'w')
  stats_file.writelines('cutoff;' + ';'.join(header))
  stats_file.write('\n')
  cutoffs = sorted(cutoff_stats.keys())

  for cutoff in cutoffs:
    values = cutoff_stats[cutoff]
    if values:
      stats_file.write('%.3f;' % cutoff)
      for key in header:
        stats_file.write('%.3f;' % values[key])
    stats_file.write('\n')

  stats_file.close()
  # Display on the console
  if args.verbose:
    print '\nCutoff stats\n'
    sys.stdout.writelines('cutoff;' + ';'.join(header))
    sys.stdout.write('\n')

    for cutoff in cutoffs:
      values = cutoff_stats[cutoff]
      if values:
        sys.stdout.write('%.3f;' % cutoff)
        for key in header:
          sys.stdout.write('%.3f;' % values[key])
      sys.stdout.write('\n')

print 'The log file saved to: %s' % log_file
print 'See you next time!'
logger.close()
