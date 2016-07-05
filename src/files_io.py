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

BondMatcher - I/O library. Handles opening and writing different files."""


from collections import defaultdict
import numpy
import os

import structures  # pylint:disable=W0403


def prepare_path(file_path):
  """Prepare the file to open.

  Args:
    file_path: The file path.

  Returns:
    The path to the file.
  """

  if os.path.exists(file_path):
    file_name = os.path.basename(file_path)
    dir_name = os.path.dirname(file_path)
    if not dir_name:
      dir_name = '.'
    existing_copies = [x for x in os.listdir(dir_name) if x.startswith('_%s' % file_name)]
    if existing_copies:
      max_copy_id = max([int(x.strip('_').split('.')[-1]) for x in existing_copies])
    else:
      max_copy_id = 0
    new_file_name = '_%s.%d_' % (file_name, max_copy_id+1)
    new_file_path = os.path.join(dir_name, new_file_name)
    print '\nFound: %s, backup on: %s\n' % (file_path, new_file_path)
    os.rename(file_path, new_file_path)

  return file_path


class File(object):
  """File object.

  Args:
    opened: The flag indicated that the file is opened.
    content: The raw content of the file.
    chains: The dict with chains.
    chain_idx_type_map: The map from idx to type.
    chain_type_idx_map: The map from type to chain idx.
    box: The tuple with simulation box dimensions.
    scale_factor: How to scale the numbers, by default we use Angstroms.
  """
  opened = False
  content = None
  chains = {}
  chain_idx_type_map = {}
  chain_type_idx_map = {}
  box = None
  data = None
  scale_factor = 1.0
  file = None


  def __init__(self, file_name):
    self.file_name = file_name

  def open(self):
    """Open the pdb file."""
    self.file = open(self.file_name, 'r')
    self.opened = True


class GROFile(File):
  scale_factor = 10.0

  def read(self):
    """Reads the .gro file and return the atom list.

    Returns:
      The dict with atoms (key: atom_id, value: atom object).
    """

    atoms = {}
    if not self.opened:
      raise Exception('File is not open.')

    if not self.content:
      self.content = self.file.readlines()

    number_of_atoms = int(self.content[1])

    for line in self.content[2:number_of_atoms + 2]:
      chain_idx = int(float(line[0:5].strip()))
      chain_name = line[5:10].strip()
      at_name = line[10:15].strip()
      at_id = int(line[15:20].strip())
      # Nedd to rescale.
      pos_x = float(line[20:28].strip()) * self.scale_factor
      pos_y = float(line[28:36].strip()) * self.scale_factor
      pos_z = float(line[36:44].strip()) * self.scale_factor

      atoms[at_id] = (
          structures.Atom(
              atom_id=at_id,
              name=at_name,
              chain_name=chain_name,
              chain_idx=chain_idx,
              position=numpy.array([pos_x, pos_y, pos_z])
          ))

      if chain_idx not in self.chains:
        self.chains[chain_idx] = set()
        self.chain_idx_type_map[chain_idx] = chain_name
      self.chains[chain_idx].add(atoms[at_id])

    # Reads the box size, the last line.
    self.box = numpy.array(
        map(float, filter(None, self.content[number_of_atoms + 2].split(' ')))
        ) * self.scale_factor

    self.data = atoms


class PDBFile(File):

  def read(self):
    """Reads the file and return atom list."""
    atoms = {}
    if not self.opened:
      raise Exception('File is not open.')

    if not self.content:
      self.content = self.file.readlines()

    for line in self.content:
      if line.startswith('CRYST1'):
        # Box size
        self.box = numpy.array(
            map(float, filter(None, line.split(' '))[1:4])
            )
      elif line.startswith('ATOM') or line.startswith('HETATM'):
        atom_id = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        chain_name = line[17:20].strip()  # Residue name
        chain_idx = line[22:26].strip()
        pos_x = float(line[30:38])
        pos_y = float(line[38:46])
        pos_z = float(line[46:54])
        atoms[atom_id] = (
          structures.Atom(
              atom_id=atom_id,
              name=atom_name,
              chain_name=chain_name,
              chain_idx=chain_idx,
              position=numpy.array([pos_x, pos_y, pos_z])
              ))
        if chain_idx not in self.chains:
          self.chains[chain_idx] = set()
          self.chain_idx_type_map[chain_idx] = chain_name
        self.chains[chain_idx].add(atoms[atom_id])

    self.data = atoms


  def write(self, output_file, line_replace=None):
    """Write the file again."""
    if not self.opened:
      raise Exception('File is not open.')

    for line in self.content:
      line_replacer = line_replace.get(line[0:6].strip())
      if line_replacer:
        output_file.write(line_replacer(line))
      else:
        output_file.write(line)


class TopologyFile(File):
  """Very basic representation of topology file."""

  # Store the new data.
  new_data = {
      'bonds': {},
      'angles': [],
      'dihedrals': [],
      'pairs': []
      }
  current_section = None
  # Some helper structures for bonds.
  # key: atom_i
  # val: set of atom_j
  bonds_def = defaultdict(set)

  # Current set of data.
  bonds = []
  atoms = {}
  chains = {}
  chain_idx = {}
  angles = []
  dihedrals = []
  pairs = []

  # Idx.
  _indexes = {
      'bonds': 0,
      'angles': 0,
      'dihedrals': 0,
      'pairs': 0
      }

  def __init__(self, file_name):
    super(TopologyFile, self).__init__(file_name)
    self.parsers = {
        'atoms': self._parse_atoms,
        'bonds': self._parse_bonds,
        'dihedrals': self._parse_dihedrals,
        'pairs': self._parse_pairs,
        'angles': self._parse_angles
        }

    self.writers = {
        'bonds': self._write_bonds,
        'angles': self._write_default,
        'dihedrals': self._write_default,
        'pairs': self._write_pairs
        }

  def add_new_bonds(self, new_bonds):
    """Adds the new bonds to the current set.

    Args:
      new_bonds: The dict with the key as atom_pair and value as the list
        of new bonds
    """
    for bonds in new_bonds.values():
      for bond in bonds:
        self.data['bonds'][bond] = [3]
        if bond[0] not in self.bonds_def:
          self.bonds_def[bond[0]] = set()
        self.bonds_def[bond[0]].add(bond[1])
        if bond[1] not in self.bonds_def:
          self.bonds_def[bond[1]] = set()
        self.bonds_def[bond[1]].add(bond[0])

  # Parsers for the data.
  def _parse_bonds(self, raw_data):
    atom_tuple = tuple(map(int, raw_data[0:2]))
    self.bonds.append(atom_tuple)
    self._indexes['bonds'] += 1

    self.bonds_def[atom_tuple[0]].add(atom_tuple[1])
    self.bonds_def[atom_tuple[1]].add(atom_tuple[0])

  def _parse_atoms(self, raw_data):
    at = structures.Atom(
        atom_id=int(raw_data[0]),
        name=raw_data[4],
        chain_name=raw_data[3],
        chain_idx=int(raw_data[2])
        )
    if at.chain_name not in self.chains:
      self.chains[at.chain_name] = defaultdict(list)

    if at.chain_name not in self.chain_idx:
      self.chain_idx[at.chain_name] = {}
    if at.chain_idx not in self.chain_idx[at.chain_name]:
      self.chain_idx[at.chain_name][at.chain_idx] = {}

    self.chains[at.chain_name][at.chain_idx].append(at)
    self.chain_idx[at.chain_name][at.chain_idx][at.atom_id] = at
    self.atoms[at.atom_id] = at

  def _parse_angles(self, raw_data):
    self.angles.append(tuple(map(int, raw_data[0:3])))
    self._indexes['angles'] += 1

  def _parse_dihedrals(self, raw_data):
    self.dihedrals.append(tuple(map(int, raw_data[0:4])))
    self._indexes['dihedrals'] += 1

  def _parse_pairs(self, raw_data):
    self.pairs.append(tuple(map(int, raw_data[0:2])))
    self._indexes['pairs'] += 1

  # Writers
  def _write_bonds(self, data):  # pylint:disable=R0201
    return_data = []
    for key, values in data.iteritems():
      if (key[0], key[1]) not in self.bonds or (key[1], key[0]) not in self.bonds:
        return_data.append('%d %d %s' % (
          key[0], key[1], ' '.join(map(str, values))
          ))
    return return_data

  def _write_pairs(self, data):  # pylint:disable=R0201
    return_data = []
    for key, values in data.iteritems():
      if (key[0], key[1]) not in self.pairs or (key[1], key[0]) not in self.pairs:
        return_data.append('%d %d %s' % (
          key[0], key[1], ' '.join(map(str, values))
          ))
    return return_data

  def _write_default(self, data):  # pylint:disable=R0201
    return ['%s' % ' '.join(map(str, x)) for x in data]

  def read(self):
    """Reads the topology file."""
    if not self.opened:
      raise Exception('File is not open.')

    if not self.content:
      self.content = self.file.readlines()

    # New version
    current_parser = None
    visited_sections = set()
    section_name = None
    previous_section = None
    for line in self.content:
      line = line.strip()
      if line.startswith(';') or line.startswith('#'):
        continue
      elif line.startswith('['):  # Section
        previous_section = section_name
        section_name = line.replace('[', '').replace(']', '').strip()
        current_parser = self.parsers.get(section_name)
        visited_sections.add(previous_section)
      else:
        if current_parser is not None and section_name not in visited_sections:
          raw_data = filter(None, line.split(' '))
          if raw_data:
              current_parser(raw_data)  # pylint:disable=E1102
	  else:
              if section_name in self._indexes: 
                  self._indexes[section_name] += 1
    self.cleanup()

  def cleanup(self):
    pass

  def write(self, filename=None):
    """Updates the topology file.

    Args:
      filename: The optional output filename.
    """

    if not self.opened:
      raise Exception('The file is not open.')

    if filename:
      output_file = open(prepare_path(filename), 'w')
    else:
      output_file = open(prepare_path(self.file_name), 'w')

    new_data = []
    current_section = None
    previous_section = None
    visited_sections = set()
    line_idx = 0
    writer = None
    for line in self.content:
      tmp_line = line.strip()
      if tmp_line.startswith('['):
        new_data.append(line)
        previous_section = current_section
        current_section = tmp_line.replace('[', '').replace(']', '').strip()
        visited_sections.add(previous_section)
        writer = self.writers.get(current_section, lambda x: x)
        line_idx = 0
      else:
        new_data.append(line)
        line_idx += 1
        if (line_idx == self._indexes.get(current_section, -1)
            and current_section not in visited_sections):
          # We have to add the new data.
          new_data.append('; new %s\n' % current_section)
          new_data.extend([
            '%s\n' % x for x in writer(self.new_data.get(current_section, []))  # pylint:disable=E1102
            ])
    output_file.writelines(new_data)
    output_file.close()
