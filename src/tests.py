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

BondMatcher - unit tests."""

import numpy
import unittest

# pylint: disable=W0403
import analysis
import files_io
import structures
import tools

# Some common settings
class TestSettings(structures.BondSettings):
  CHAIN_1 = 'ch1'
  CHAIN_2 = 'ch2'
  # Define some constants, name, chain_name, max_degree

  CH1_Z1 = structures.Atom(name='Z1', chain_name=CHAIN_1, max_degree=3)
  CH1_Z2 = structures.Atom(name='Z2', chain_name=CHAIN_1, max_degree=3)
  CH1_Z3 = structures.Atom(name='Z3', chain_name=CHAIN_1, max_degree=3)

  CH2_C1 = structures.Atom(name='C1', chain_name=CHAIN_2, max_degree=3)
  CH2_C2 = structures.Atom(name='C2', chain_name=CHAIN_2, max_degree=3)
  CH2_H1 = structures.Atom(name='H1', chain_name=CHAIN_2, max_degree=1)

  BOND_1 = (CH1_Z1, CH2_C1)

  # The atom symbols that should form a bond.
  ATOM_PAIRS = [BOND_1]

  # For given atom pairs define the length of the bond and the frequency (according to the
  # selected force field). Please use the same order as in ATOM_PAIRS
  BOND_CONFIG = {
    BOND_1: [1, 0.146, 319657.6],
    }

  # Setup for the angles, order does matter - it's a pattern not set.
  ANGLES = {
    BOND_1: {
        (CH2_C2, CH2_C1, CH1_Z1): [],
        (CH1_Z2, CH1_Z1, CH2_C1): [],
        (CH1_Z1, CH2_C1, CH1_Z1): []
        }
    }

  DIHEDRALS = {
    BOND_1: {
        (CH2_C2, CH2_C1, CH1_Z1, CH1_Z2): [],
        (CH2_H1, CH2_C2, CH2_C1, CH1_Z1): [],
        (CH1_Z3, CH1_Z2, CH1_Z1, CH2_C1): [],
        (CH1_Z1, CH2_C1, CH1_Z1, CH1_Z2): []
        }
    }
  CUTOFF = 1


class BaseTestCase(unittest.TestCase):

  def assertSameElements(self, one, two):
    """Checks if two iterable elements are the same. The order does not matter."""
    missing = []

    for el_one in one:
      if el_one not in two:
        missing.append(el_one)

    self.assertEqual(missing, [])


class TestBaseTestCase(BaseTestCase):

  def testAssertSameElementsOk(self):
    one = [1, 2, 3, 4]
    two = [3, 4, 1, 2]

    self.assertSameElements(one, two)

  def testAsssertsameElementsFail(self):
    one = [1, 2, 3]
    two = [1, 2, 4]

    try:
      self.assertSameElements(one, two)
    except Exception as ex:
      self.assertIsInstance(ex, AssertionError)



class TestAnalysis(BaseTestCase):

  def setUp(self):
    super(TestAnalysis, self).setUp()
    self.settings = TestSettings()

    self.all_atoms = {}
    # Build atoms.
    chains = {
        'ch1': [
            ((1, 'C1', 1, 0), (2, 'Z3', 1, 1), (3, 'C2', 0, 1), (4, 'Z2', 1, 2), (5, 'Z1', 2, 2)),
            ((6, 'C1', 5, 6), (7, 'Z3', 5, 5), (8, 'C2', 4, 5), (9, 'Z2', 5, 4), (10, 'Z1', 5, 3))
            ],
        'ch2': [((11, 'C1', 3, 2), (12, 'C2', 4, 2), (13, 'H1', 4, 1), (14, 'C1', 5, 2))]
        }
    bonds = [
        (1, 2), (2, 3), (2, 4), (4, 5),
        (6, 7), (7, 8), (7, 9), (9, 10),
        (11, 12), (12, 13), (12, 14)
        ]
    for chain_name, chs in chains.iteritems():
      for ch in chs:
        for at in ch:
          self.all_atoms[at[0]] = structures.Atom(
              atom_id=at[0],
              name=at[1],
              chain_name=chain_name,
              position=numpy.array([at[2], at[3], 0])
              )

    self.top_file = files_io.TopologyFile('dummy')
    self.top_file.new_data['bonds'] = {}
    self.top_file.bonds_def = {}
    for bond in bonds:
      self.top_file.new_data['bonds'][bond] = [3]
      if bond[0] not in self.top_file.bonds_def:
        self.top_file.bonds_def[bond[0]] = set()
      self.top_file.bonds_def[bond[0]].add(bond[1])
      if bond[1] not in self.top_file.bonds_def:
        self.top_file.bonds_def[bond[1]] = set()
      self.top_file.bonds_def[bond[1]].add(bond[0])

    self.valid_atoms = None
    self.number_of_bonds = 0
    self.valid_atoms, self.number_of_bonds = tools.prepare_atoms(
        self.all_atoms, self.settings, self.top_file)

  def test_find_atom_pairs(self):
    new_bonds = analysis.find_atom_pairs(
        self.valid_atoms, self.settings.CUTOFF, self.settings.ATOM_PAIRS)
    self.assertEqual(len(new_bonds.values()[0]), 2)
    new_bond_list = new_bonds[self.settings.ATOM_PAIRS[0]]
    self.assertEqual(
        [(10, 14), (11, 5)],
        [(x[0].atom_id, x[1].atom_id) for x in new_bond_list])
    # Check the degree of the atoms, should increase
    self.assertEqual([x.degree for x in self.valid_atoms], [2, 2, 2, 2])

  def test_angles(self):
    new_bonds = analysis.find_atom_pairs(
        self.valid_atoms, self.settings.CUTOFF, self.settings.ATOM_PAIRS)
    self.top_file.add_new_bonds(new_bonds)
    find_angles = analysis.Analysis(
        new_bonds, self.settings.ANGLES, self.settings.DIHEDRALS)

    new_angles = find_angles.find_angles()
    self.assertSameElements(
        new_angles,
        [[9, 10, 14], [4, 5, 11], [12, 14, 10], [12, 11, 5]]
        )
    new_dihedrals = find_angles.find_dihedrals()
    dih = [[13, 12, 14, 10],
          [12, 14, 10, 9],
          [2, 4, 5, 11],
          [7, 9, 10, 14],
          [13, 12, 11, 5],
          [12, 11, 5, 4]]
    self.assertEqual(len(dih), len(new_dihedrals))
    self.assertListEqual(new_dihedrals, dih)


class TestAnalysisTwoChain(BaseTestCase):

  def setUp(self):
    super(TestAnalysisTwoChain, self).setUp()
    self.settings = TestSettings()

    self.all_atoms = {}
    # Build atoms.
    chains = {
        'ch1': [
            ((1, 'C1', 1, 0), (2, 'Z3', 1, 1), (3, 'C2', 0, 1), (4, 'Z2', 1, 2), (5, 'Z1', 2, 2)),
            ((6, 'C1', 5, 6), (7, 'Z3', 5, 5), (8, 'C2', 4, 5), (9, 'Z2', 5, 4), (10, 'Z1', 5, 3)),
            ((15, 'C1', 8, 1), (16, 'Z3', 7, 1), (17, 'C2', 7, 2),
             (18, 'Z2', 6, 1), (19, 'Z1', 5, 1)),
            ],
        'ch2': [((11, 'C1', 3, 2), (12, 'C2', 4, 2), (13, 'H1', 4, 1), (14, 'C1', 5, 2))]
        }
    bonds = [
        (1, 2), (2, 3), (2, 4), (4, 5),
        (6, 7), (7, 8), (7, 9), (9, 10),
        (11, 12), (12, 13), (12, 14),
        (15, 16), (16, 17), (16, 18), (18, 19)
        ]
    for chain_name, chs in chains.iteritems():
      for ch in chs:
        for at in ch:
          self.all_atoms[at[0]] = structures.Atom(
              atom_id=at[0],
              name=at[1],
              chain_name=chain_name,
              position=numpy.array([at[2], at[3], 0])
              )

    self.top_file = files_io.TopologyFile('dummy')
    self.top_file.new_data['bonds'] = {}
    self.top_file.bonds_def = {}
    for bond in bonds:
      self.top_file.new_data['bonds'][bond] = [3]
      if bond[0] not in self.top_file.bonds_def:
        self.top_file.bonds_def[bond[0]] = set()
      self.top_file.bonds_def[bond[0]].add(bond[1])
      if bond[1] not in self.top_file.bonds_def:
        self.top_file.bonds_def[bond[1]] = set()
      self.top_file.bonds_def[bond[1]].add(bond[0])

    self.valid_atoms = None
    self.number_of_bonds = 0
    self.valid_atoms, self.number_of_bonds = tools.prepare_atoms(
        self.all_atoms, self.settings, self.top_file)

  def test_find_atom_pairs(self):
    new_bonds = analysis.find_atom_pairs(
        self.valid_atoms, self.settings.CUTOFF, self.settings.ATOM_PAIRS)
    self.assertEqual(len(new_bonds.values()[0]), 3)
    new_bond_list = new_bonds[self.settings.ATOM_PAIRS[0]]
    self.assertListEqual(
        [(19, 14), (10, 14), (11, 5)],
        [(x[0].atom_id, x[1].atom_id) for x in new_bond_list])
    # Check the degree of the atoms, should increase
    self.assertEqual([x.degree for x in self.valid_atoms], [2, 2, 2, 2, 3])

  def test_angles(self):
    new_bonds = analysis.find_atom_pairs(
        self.valid_atoms, self.settings.CUTOFF, self.settings.ATOM_PAIRS)
    self.top_file.add_new_bonds(new_bonds)
    find_angles = analysis.Analysis(
        new_bonds, self.settings.ANGLES, self.settings.DIHEDRALS)

    new_angles = find_angles.find_angles()
    self.assertSameElements(
        new_angles,
        [[9, 10, 14], [4, 5, 11], [12, 14, 10], [12, 11, 5],
         [18, 19, 14], [10, 14, 19], [12, 14, 19]]
        )

  def test_dihedrals(self):
    new_bonds = analysis.find_atom_pairs(
        self.valid_atoms, self.settings.CUTOFF, self.settings.ATOM_PAIRS)
    self.top_file.add_new_bonds(new_bonds)
    find_angles = analysis.Analysis(
        new_bonds, self.settings.ANGLES, self.settings.DIHEDRALS)

    new_dihedrals = find_angles.find_dihedrals()
    self.assertSameElements(
        new_dihedrals,
        [[7, 9, 10, 14],
         [2, 4, 5, 11],
         [13, 12, 14, 10],
         [13, 12, 11, 5],
         [12, 14, 10, 9],
         [12, 11, 5, 4],
         [12, 14, 19, 18],
         [13, 12, 14, 19],
         [16, 18, 19, 14],
         [19, 14, 10, 9]
        ]
        )


if __name__ == '__main__':
  unittest.main()
