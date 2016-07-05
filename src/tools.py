"""
Copyright (C) 2014
  KU Leuven, Computer Science Department
  Jakub Krajniak <jkrajniak@gmail.com>

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

"""

from collections import defaultdict
import datetime
import os
import sys


class Logger(object):
  def __init__(self, root_dir, log_name='run_log.log'):
    self.terminal = sys.stdout
    self.log_path = os.path.join(root_dir, log_name)
    self.log = open(self.log_path, 'w')
    self.log.write('LOG FILE date=%s\n' % datetime.datetime.now())
    self.log.write('PWD: %s\n' % os.environ.get('PWD'))
    self.log.write('USER: %s\n============================\n' % os.environ.get('USER'))

  def close(self):
    self.log.write(
        ('============================\n'
         'LOG file closed date=%si\n') % datetime.datetime.now())
    self.log.close()

  def write(self, message):
    self.terminal.write(message)
    self.log.write(message)

  def writelines(self, lines):
    for l in lines:
      self.terminal.write(l)
      self.log.write(l)

  def write_to_file(self, message):
    self.log.write(message)

  def enable(self):
    sys.stdout = self

  def disable(self):
    sys.stdout = self.terminal


def update_degree(atoms, bonds):
  """Update the degree values for all of atoms.

  Args:
    atoms: The list of all atoms of the system.
    bonds: The dictionary with the atom_i: atom_j structure

  Returns:
    The updates list of all atoms.
  """
  bond_ids = [x[0] for x in bonds]
  bond_ids.extend([x[1] for x in bonds])

  for atom in atoms:
    atom.degree += bond_ids.count(atom.atom_id)

  return atoms


def find_neighbours(atoms, top_file, all_atoms):
  """Find the neighbours of the atoms.

  Args:
    atoms: The list of atoms to update.
    top_file: The topology file.
    all_atoms: The dict with all atoms in the system.

  Returns:
    The update list of atoms.
  """
  for at in atoms:
    at.neighbours = [
        all_atoms[j] for j in top_file.bonds_def.get(at.atom_id, [])
        ]

  return atoms


def get_valid_atoms_by_name(atoms, atom_pairs):
  """Returns only atoms that has specific name.

  Args:
    atoms: The dict with all atoms in the system.
    atom_pairs: The description of valid atom symbols, it's a list of tuples.

  Returns:
    The list of atoms.
  """
  # Help structure.
  chain_ids_atoms = {}
  for pairs in atom_pairs:
    for p in pairs:
      if p.chain_name not in chain_ids_atoms:
        chain_ids_atoms[p.chain_name] = []
      chain_ids_atoms[p.chain_name].append(p)

  return filter(
      lambda x: x in chain_ids_atoms.get(x.chain_name, []),
      atoms.values()
      )


def get_valid_atoms_by_bonds(atoms, atom_pairs, bonds):
  """Returns only atoms that not formed bonds.

  Args:
    atoms: The list of atoms.
    atom_pairs: The description of bonds.
    bonds: The current set of bonds.

  Returns:
    The tuple with list of atoms, the number of bonds already formed
  """
  # Help structure.
  atoms_def_per_chain = {}
  valid_bonds = {}
  for pairs in atom_pairs:
    valid_bonds[pairs] = []
    for p in pairs:
      if p.chain_name not in atoms_def_per_chain:
        atoms_def_per_chain[p.chain_name] = {}
      atoms_def_per_chain[p.chain_name][p.name] = p

  valid_atoms = {x.atom_id: x for x in atoms}
  number_of_bonds = 0

  for atom_id_i, atom_id_j in bonds:
    atom_i = valid_atoms.get(atom_id_i)
    atom_j = valid_atoms.get(atom_id_j)

    if atom_i and atom_j:
      key1 = (atom_i, atom_j)
      key2 = (atom_j, atom_i)
     
      if key1 in valid_bonds:
        number_of_bonds += 1
        valid_bonds[key1].append(key1)
      elif key2 in valid_bonds:
        number_of_bonds += 1
        valid_bonds[key2].append(key1)

  return valid_atoms.values(), number_of_bonds, valid_bonds


def get_cutoff(density, cutoff_settings, default):
  """Returns the cutoff, based on the density.

  Args:
    density: The number.
    cutoff_settings: The dict with keys as the minimum density and the value as the cutoff.
    default: The default cutoff.

  Returns:
    The cutoff distance.
  """
  return cutoff_settings.get(filter(lambda x: x <= density, sorted(cutoff_settings))[-1], default)


def find_chains(top_file, pdb_file, settings):
  """Find the chains that are involved in the bonds.

  Args:
    top_file: The top file.
    pdb_file: The pdb file.
    settings: The settings

  Returns:
    The tuple with chain_ids involved in the network and chain_ids involved in the bonds, divided
    by the chain name.
  """
  chain_types_in_bonds = defaultdict(set)
  chain_ids_in_the_bonds = set()

  chain_pairs = set()
  for b in settings.ATOM_PAIRS:
    chain_pairs.add((b[0].chain_name, b[1].chain_name))

  # Lets check the bonds
  for bond in top_file.bonds:
    ai = pdb_file.data[bond[0]]
    aj = pdb_file.data[bond[1]]
    chain_pair_1 = (ai.chain_name, aj.chain_name)
    chain_pair_2 = (aj.chain_name, ai.chain_name)
    if chain_pair_1 in chain_pairs or chain_pair_2 in chain_pairs:
      chain_ids_in_the_bonds.update([ai.chain_idx, aj.chain_idx])
      chain_types_in_bonds[ai.chain_name].add(aj.chain_idx)
      chain_types_in_bonds[aj.chain_name].add(ai.chain_idx)

  return chain_ids_in_the_bonds, chain_types_in_bonds


def prepare_atoms(all_atoms, settings, top_file):
  """Run the chain of functions to prepare the atoms to process.

  Args:
    all_atoms: The dict with all atoms.
    settings: The settings object.
    top_file: The top_file.

  Returns:
    The tuple with list of valid atoms and number of bonds.
  """

  important_atoms = [x for x in all_atoms.values() if x in settings.IMPORTANT_ATOMS]
  # Attache the neighbours to the atoms.
  find_neighbours(important_atoms, top_file, all_atoms)

  valid_atoms = get_valid_atoms_by_name(all_atoms, settings.ATOM_PAIRS)
  valid_atoms = update_degree(valid_atoms, top_file.bonds)
  valid_atoms, number_of_bonds, valid_bonds = get_valid_atoms_by_bonds(
      valid_atoms,
      settings.ATOM_PAIRS,
      top_file.bonds
      )

  for v in important_atoms:
    k = (v.name, v.chain_name)
    at_def = settings.ATOMS_DEF[k]
    v.max_degree = at_def.max_degree
    v.radius = at_def.radius

  return valid_atoms, number_of_bonds, valid_bonds
