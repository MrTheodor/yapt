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

BondMatcher - unit tests for the analysis part
"""

import itertools
import unittest

# pylint: disable=W0403
import analysis
import structures
import test_case


# Some common settings
class TestSettings(structures.BondSettings):
  CHAIN_1 = 'ch1'
  CHAIN_2 = 'ch2'
  # Define some constants, name, chain_name, max_degree

  CH1_Z1 = structures.Atom(name='Z1', chain_name=CHAIN_1, max_degree=3)
  CH1_Z2 = structures.Atom(name='Z2', chain_name=CHAIN_1, max_degree=3)
  CH1_Z3 = structures.Atom(name='Z3', chain_name=CHAIN_1, max_degree=3)
  CH1_C1 = structures.Atom(name='C1', chain_name=CHAIN_1, max_degree=4)
  CH1_C2 = structures.Atom(name='C2', chain_name=CHAIN_1, max_degree=4)
  CH1_H11 = structures.Atom(name='H11', chain_name=CHAIN_1, max_degree=1)
  CH1_H12 = structures.Atom(name='H12', chain_name=CHAIN_1, max_degree=1)
  CH1_H21 = structures.Atom(name='H21', chain_name=CHAIN_1, max_degree=1)
  CH1_H22 = structures.Atom(name='H22', chain_name=CHAIN_1, max_degree=1)

  CH2_C1 = structures.Atom(name='C1', chain_name=CHAIN_2, max_degree=3)
  CH2_C2 = structures.Atom(name='C2', chain_name=CHAIN_2, max_degree=3)
  CH2_C3 = structures.Atom(name='C3', chain_name=CHAIN_2, max_degree=3)
  CH2_H11 = structures.Atom(name='H11', chain_name=CHAIN_2, max_degree=1)
  CH2_H12 = structures.Atom(name='H12', chain_name=CHAIN_2, max_degree=1)
  CH2_O1 = structures.Atom(name='O1', chain_name=CHAIN_2, max_degree=2)
  CH2_H2 = structures.Atom(name='H2', chain_name=CHAIN_2, max_degree=2)

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
        (CH2_H11, CH2_C2, CH2_C1, CH1_Z1): [],
        (CH1_Z3, CH1_Z2, CH1_Z1, CH2_C1): [],
        (CH1_Z1, CH2_C1, CH1_Z1, CH1_Z2): []
        }
    }

  PAIRS = {
      BOND_1: [1]
      }


class TestFind14Pairs(test_case.BaseTestCase):

  def setUp(self):
    super(TestFind14Pairs, self).setUp()
    self.conf = TestSettings()
    self.new_bonds = {}
    self.all_atoms = {}

  def _create_simple(self):
    conf = self.conf
    atom_id = itertools.count(1)
    # Creates linker chain with both active two ends.
    create_atom = structures.Atom.create
    ch1z11 = structures.Atom.create(conf.CH1_Z1, next(atom_id), 1)
    ch1c11 = structures.Atom.create(conf.CH1_C1, next(atom_id), 1)
    ch1c21 = structures.Atom.create(conf.CH1_C2, next(atom_id), 1)
    ch1h11 = structures.Atom.create(conf.CH1_H11, next(atom_id), 1)
    ch1h12 = structures.Atom.create(conf.CH1_H12, next(atom_id), 1)

    ch1z12 = structures.Atom.create(conf.CH1_Z1, next(atom_id), 1)
    ch1c22 = structures.Atom.create(conf.CH1_C2, next(atom_id), 1)
    ch1c12 = structures.Atom.create(conf.CH1_C1, next(atom_id), 1)
    ch1h21 = structures.Atom.create(conf.CH1_H21, next(atom_id), 1)
    ch1h22 = structures.Atom.create(conf.CH1_H22, next(atom_id), 1)
    # Create the nb lists.
    ch1z11.neighbours = [ch1c11]
    ch1c11.neighbours = [ch1z11, ch1h11, ch1h12, ch1c21]
    ch1h12.neighbours = [ch1c11]
    ch1h11.neighbours = [ch1c11]
    ch1h21.neighbours = [ch1c11]
    ch1z12.neighbours = [ch1c22]
    ch1c22.neighbours = [ch1z12, ch1h21, ch1h22, ch1c12]
    ch1h21.neighbours = [ch1c22]
    ch1h22.neighbours = [ch1c22]
    ch1c12.neighbours = [ch1c22]
    all_atoms = [
        ch1z11, ch1c11, ch1c21, ch1h11, ch1h12, ch1z12, ch1c22, ch1c12, ch1h21, ch1h22
        ]
    bonds = []
    # Create the polymer chains.
    for chain_id in [2]:
      c1 = create_atom(conf.CH2_C1, next(atom_id), chain_id)
      c2 = create_atom(conf.CH2_C2, next(atom_id), chain_id)
      c3 = create_atom(conf.CH2_C3, next(atom_id), chain_id)
      h11 = create_atom(conf.CH2_H11, next(atom_id), chain_id)
      h12 = create_atom(conf.CH2_H12, next(atom_id), chain_id)
      h2 = create_atom(conf.CH2_H2, next(atom_id), chain_id)
      o1 = create_atom(conf.CH2_O1, next(atom_id), chain_id)
      # Update the nb lists.
      c1.neighbours = [h11, c2, h12]
      c2.neighbours = [c1, h2, o1, c3]
      c3.neighbours = [c2]
      h12.neighbours = [c1]
      h11.neighbours = [c1]
      h2.neighbours = [c2]
      o1.neighbours = [c2]
      all_atoms.extend([
        c1, c2, c3, h11, h12, h2, o1
        ])
      if chain_id in (2, 3):
        ch1z11.neighbours.append(c1)
        c1.neighbours.append(ch1z11)
        bonds.append((c1, ch1z11))
      else:
        ch1z12.neighbours.append(c1)
        c1.neighbours.append(ch1z12)
        bonds.append((c1, ch1z12))

    self.all_atoms = {x.atom_id: x for x in all_atoms}
    self.new_bonds = {conf.BOND_1: bonds}

  def _create_complex(self):
    conf = self.conf
    atom_id = itertools.count(1)
    # Creates linker chain with both active two ends.
    create_atom = structures.Atom.create
    ch1z11 = structures.Atom.create(conf.CH1_Z1, next(atom_id), 1)
    ch1c11 = structures.Atom.create(conf.CH1_C1, next(atom_id), 1)
    ch1c21 = structures.Atom.create(conf.CH1_C2, next(atom_id), 1)
    ch1h11 = structures.Atom.create(conf.CH1_H11, next(atom_id), 1)
    ch1h12 = structures.Atom.create(conf.CH1_H12, next(atom_id), 1)

    ch1z12 = structures.Atom.create(conf.CH1_Z1, next(atom_id), 1)
    ch1c22 = structures.Atom.create(conf.CH1_C2, next(atom_id), 1)
    ch1c12 = structures.Atom.create(conf.CH1_C1, next(atom_id), 1)
    ch1h21 = structures.Atom.create(conf.CH1_H21, next(atom_id), 1)
    ch1h22 = structures.Atom.create(conf.CH1_H22, next(atom_id), 1)
    # Create the nb lists.
    ch1z11.neighbours = [ch1c11]
    ch1c11.neighbours = [ch1z11, ch1h11, ch1h12, ch1c21]
    ch1h12.neighbours = [ch1c11]
    ch1h11.neighbours = [ch1c11]
    ch1h21.neighbours = [ch1c11]
    ch1z12.neighbours = [ch1c22]
    ch1c22.neighbours = [ch1z12, ch1h21, ch1h22, ch1c12]
    ch1h21.neighbours = [ch1c22]
    ch1h22.neighbours = [ch1c22]
    ch1c12.neighbours = [ch1c22]
    all_atoms = [
        ch1z11, ch1c11, ch1c21, ch1h11, ch1h12, ch1z12, ch1c22, ch1c12, ch1h21, ch1h22
        ]
    bonds = []
    # Create the polymer chains.
    for chain_id in [2, 3, 4]:
      c1 = create_atom(conf.CH2_C1, next(atom_id), chain_id)
      c2 = create_atom(conf.CH2_C2, next(atom_id), chain_id)
      c3 = create_atom(conf.CH2_C3, next(atom_id), chain_id)
      h11 = create_atom(conf.CH2_H11, next(atom_id), chain_id)
      h12 = create_atom(conf.CH2_H12, next(atom_id), chain_id)
      h2 = create_atom(conf.CH2_H2, next(atom_id), chain_id)
      o1 = create_atom(conf.CH2_O1, next(atom_id), chain_id)
      # Update the nb lists.
      c1.neighbours = [h11, c2, h12]
      c2.neighbours = [c1, h2, o1, c3]
      c3.neighbours = [c2]
      h12.neighbours = [c1]
      h11.neighbours = [c1]
      h2.neighbours = [c2]
      o1.neighbours = [c2]
      all_atoms.extend([
        c1, c2, c3, h11, h12, h2, o1
        ])
      if chain_id in (2, 3):
        ch1z11.neighbours.append(c1)
        c1.neighbours.append(ch1z11)
        bonds.append((c1, ch1z11))
      else:
        ch1z12.neighbours.append(c1)
        c1.neighbours.append(ch1z12)
        bonds.append((c1, ch1z12))

    self.new_bonds = {conf.BOND_1: bonds}
    self.all_atoms = {x.atom_id: x for x in all_atoms}

  def testBasic(self):
    self._create_simple()
    new_pairs = analysis.find_1_4_pairs(self.new_bonds, self.conf).keys()
    valid_pairs = [
        [1, 13],
        [1, 16],
        [1, 17],
        [2, 12],
        [2, 14],
        [2, 15],
        [3, 11],
        [4, 11],
        [5, 11]]
    self.assertEqual(len(new_pairs), len(valid_pairs))
    self.assertSameElementsNestedList(new_pairs, valid_pairs)

  def testComplex(self):
    self._create_complex()
    new_pairs = analysis.find_1_4_pairs(self.new_bonds, self.conf).keys()
    valid_pairs = [
        [1, 13],
        [1, 16],
        [1, 17],
        [2, 12],
        [2, 14],
        [2, 15],
        [3, 11],
        [4, 11],
        [5, 11],

        [1, 20],
        [1, 23],
        [1, 24],
        [2, 21],
        [2, 22],
        [2, 19],
        [3, 18],
        [4, 18],
        [5, 18],
        [19, 11],
        [18, 12],
        [18, 14],
        [18, 15],
        [21, 11],
        [22, 11],

        [6, 27],
        [6, 30],
        [6, 31],
        [7, 26],
        [7, 28],
        [7, 29],
        [8, 25],
        [9, 25],
        [10, 25]
        ]
    self.assertEqual(len(new_pairs), len(valid_pairs))
    self.assertSameElementsNestedList(new_pairs, valid_pairs)


if __name__ == '__main__':
  unittest.main()
