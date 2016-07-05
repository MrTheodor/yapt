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

BondMatcher - analysis the topology and the connectivity.
"""


from collections import defaultdict
import itertools
import numbers
import numpy
import warnings


class Analysis(object):
  """Tools for finding the angles."""

  def __init__(self, bonds, local_settings, chains):
    """Constructor of the class.

    Args:
      bonds: The dict with the new bonds.
      local_settings: The local settings.
      chains: The dict with chains.
    """
    self.bonds = bonds
    self.all_bonds = bonds
    self.angles_settings = local_settings.ANGLES
    self.dihedrals_settings = local_settings.DIHEDRALS
    self.pairs_settings = local_settings.PAIRS
    self.box = local_settings.BOX
    self.half_box = self.box / 2.0
    self.atom_pairs = local_settings.ATOM_PAIRS
    self.chains = chains
    self.local_settings = local_settings

    # Precalculate the settings
    self.atoms_def_per_chain = {}
    self.atom_pairs_settings = {}
    for atom_pair in self.atom_pairs:
      for p in atom_pair:
        if p.chain_name not in self.atoms_def_per_chain:
          self.atoms_def_per_chain[p.chain_name] = {}
        self.atoms_def_per_chain[p.chain_name][p.name] = p
      self.atom_pairs_settings[
          tuple([
              atom_pair[0].name, atom_pair[0].chain_name,
              atom_pair[1].name, atom_pair[1].chain_name
              ])] = atom_pair
      self.atom_pairs_settings[tuple([
              atom_pair[1].name, atom_pair[1].chain_name,
              atom_pair[0].name, atom_pair[0].chain_name
              ])] = atom_pair

  def add_bonds(self, bonds):
    if not self.all_bonds:
      self.all_bonds = bonds
    else:
      for k in bonds:
        self.all_bonds[k].extend(bonds[k])

  def find_angles(self):
    """Prepares angles definition.

    Returns:
      The tuples with angles definition.
    """
    new_angles = set()

    for atom_pair, atoms in self.bonds.iteritems():
      angle_config = self.angles_settings[atom_pair]

      # For every new bond create set of angles.
      for ai, aj in atoms:
        for angle_def, angle_properties in angle_config.iteritems():
          # Check where is a central atom of the bond.
          # TODO: simplify this part, looks awful!
          if ai == angle_def[1]:
            if aj == angle_def[0]:  # (aj, ai, ak)
              # Looking for the ak in the list of ai neighbours, exlcude the aj atom.
              for ak in ai.neighbours:
                if ak == angle_def[2] and ak.atom_id != aj.atom_id:
                  new_angles.add(tuple([aj.atom_id, ai.atom_id, ak.atom_id] + angle_properties))
                  break
            elif aj == angle_def[2]:  # (ak, ai, aj)
              for ak in ai.neighbours:
                if ak == angle_def[0] and ak.atom_id != aj.atom_id:
                  new_angles.add(tuple([ak.atom_id, ai.atom_id, aj.atom_id] + angle_properties))
                  break
          elif aj == angle_def[1]:
            if ai == angle_def[0]:  # (ai, aj, ak)
              for ak in aj.neighbours:
                if ak == angle_def[2] and ak.atom_id != ai.atom_id:
                  new_angles.add(tuple([ai.atom_id, aj.atom_id, ak.atom_id] + angle_properties))
                  break
            elif ai == angle_def[2]:  # (ak, aj, ai)
              for ak in aj.neighbours:
                if ak == angle_def[0] and ak.atom_id != ai.atom_id:
                  new_angles.add(tuple([ak.atom_id, aj.atom_id, ai.atom_id] + angle_properties))
                  break
    return map(list, new_angles)


  def find_dihedrals(self):
    """Prepares angles definition.

    Returns:
      The tuples with angles definition.
    """
    new_dihedrals = set()

    for atom_pair, atoms in self.bonds.iteritems():
      angle_config = self.dihedrals_settings[atom_pair]
      for at_tuple in atoms:
        for angle_def, angle_properties in angle_config.iteritems():
          at_tuple_pos = -1
          at_tuple_set = set([(x.name, x.chain_name) for x in at_tuple])
          dih = [None] * 4
          for pos in range(3):
            angle_def_set = set([(x.name, x.chain_name) for x in angle_def[pos:pos + 2]])
            if angle_def_set == at_tuple_set:
              dih[angle_def.index(at_tuple[0])] = at_tuple[0]
              dih[angle_def.index(at_tuple[1])] = at_tuple[1]
              at_tuple_pos = pos
              break
          if at_tuple_pos != -1:
            if at_tuple_pos == 0:  # 2, 3 are missing
              for nb in dih[1].neighbours:
                if nb == angle_def[2] and nb.atom_id != dih[0].atom_id:
                  dih[2] = nb
                  break
              if dih[2]:
                for nb in dih[2].neighbours:
                  if nb == angle_def[3] and nb.atom_id != dih[1].atom_id:
                    dih[3] = nb
                    break
            elif at_tuple_pos == 1:
              for nb in dih[1].neighbours:
                if nb == angle_def[0] and nb.atom_id != dih[2].atom_id:
                  dih[0] = nb
                  break
              if dih[2]:
                for nb in dih[2].neighbours:
                  if nb == angle_def[3] and nb.atom_id != dih[1].atom_id:
                    dih[3] = nb
                    break
            elif at_tuple_pos == 2:
              for nb in dih[2].neighbours:
                if nb == angle_def[1] and nb.atom_id != dih[3].atom_id:
                  dih[1] = nb
                  break
              if dih[1]:
                for nb in dih[1].neighbours:
                  if nb == angle_def[0] and nb.atom_id != dih[2].atom_id:
                    dih[0] = nb
                    break

          # TODO: simplify it, looks awful
          # Construct p1:ai p2:aj p3:ak p4:al
          if None not in dih:
            new_dihedrals.add(tuple([x.atom_id for x in dih] + angle_properties))

    return list(map(list, new_dihedrals))

  def calculate_distance(self, pos_1, pos_2):
    """Calculate distance. Minimum image convention."""
    d = pos_1 - pos_2
    for ex in range(3):
      if d[ex] > self.half_box[ex]:
        d -= self.box[ex]
      elif d[ex] < -self.half_box[ex]:
        d += self.box[ex]
    return d.dot(d)


  # Compute the angle distribution
  def unit_vector(self, vector):  # pylint:disable=R0201
    """Returns the unit vector.

    Args:
      vector: The input vector.

    Returns:
      Returns the normalized vector.
    """
    return vector / numpy.linalg.norm(vector)


  def angle_between(self, v1, v2):
    """Return the angle in radians between vectors.

    Args:
      v1: The vector a. The reference one.
      v2: The vector b.

    Returns:
      The angle between vector v1 and v2 in radians.
    """
    v1_unit = self.unit_vector(v1)
    v2_unit = self.unit_vector(v2)
    angle = numpy.arccos(numpy.dot(v1_unit, v2_unit))
    if numpy.isnan(angle):
      if (v1_unit == v2_unit).all():
        return 0.0
      else:
        return numpy.pi
    return angle


  def vector_constraint(self, atom_pair, bond):
    """Checks the angle between bond vector and the constraint vector.

    Args:
      atom_pair: The definition of the bond to use for the constraint.
      bond: The new bond that is tested.

    Returns:
      True if the new bond is valid.
    """
    confs = self.local_settings.GEOMETRIC_CONSTRAINTS['vector'].get(atom_pair, [])
    for conf in confs:
      # Localize the points, based on the definition.
      v1_begin, v1_end, v2_begin, v2_end = None, None, None, None
      v1_def = conf[0:2]
      v2_def = conf[2:4]
      if bond[0] in v1_def:
        v1_begin = bond[0]
        v2_begin = bond[1]
        v1_end = v1_def[(v1_def.index(bond[0]) - 1) % 2]
      elif bond[0] in v2_def:
        v2_begin = bond[0]
        v1_begin = bond[1]
      else:
        return None

      v1_end = v1_def[(v1_def.index(v1_begin) - 1) % 2]
      v2_end = v2_def[(v2_def.index(v2_begin) - 1) % 2]
      # Looking for the end points
      for nb in v1_begin.neighbours:
        if nb == v1_end:
          v1_end = nb
          break
      for nb in v2_begin.neighbours:
        if nb == v2_end:
          v2_end = nb
          break

      # If the vectors are not found, return None.
      if None in [
          v1_begin.position, v1_end.position, v2_begin.position, v2_end.position
          ]:
        return None

      v1 = v1_end.position - v1_begin.position
      v2 = v2_end.position - v2_begin.position
      angle = self.angle_between(v1, v2)

      if isinstance(conf[4], numbers.Number):
        return angle == numpy.radians(conf[4])
      else:
        c1 = numpy.radians(conf[4][0])
        c2 = numpy.radians(conf[4][1])
        return angle >= c1 and angle <= c2


  def sphere_constraint(self, atom_pair, bond):
    """Checks if the bond vector crosses excluded volume sphere."""
    sps = self.local_settings.GEOMETRIC_CONSTRAINTS['sphere'].get(atom_pair, [])
    spheres = []
    vb = bond[1].position - bond[0].position
    # Let's get the spheres from two molecules
    for bond_point in bond:
      for sp in sps:
        spheres.extend([x for x in self.chains[bond_point.chain_idx] if x == sp])

    # Lets analysis the spheres
    # Distance the center of the sphere to the vector should
    # be of course less than the radius.
    # http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    for sp in spheres:
      if sp.radius is None:
        raise ValueError('Please set radius for atom type: %s_%s.' % (sp.name, sp.chain_name))
      p = sp.position
      a1 = numpy.cross(bond[1].position - bond[0].position, bond[0].position - p)
      distance = a1.dot(a1)/vb.dot(vb)
      if distance <= sp.radius**2:
        return False

    return True


  def align(self, v1, v2):  # pylint:disable=R0201
    """Check if two vectors align in the same direction or not."""
    return numpy.sign(v1.dot(v2)) == 1.0


  def align_constraint(self, atom_pair, bond):  # pylint: disable=W0613,R0201
    """Checks the alignment between vectors."""
    warnings.warn('Align constraint will be implemented in the near future.')
    return True


  def cone_constraint(self, atom_pair, bond):
    """Checks if the bond vector is inside the cone.

    Args:
      atom_pair: The definition of the bond to use for the cone.
      bond: The new bond that is tested.

    Returns:
      True if the new bond is valid.
    """
    confs = self.local_settings.GEOMETRIC_CONSTRAINTS['cone'].get(atom_pair, [])
    # We need real atom points, not those definition ;) and
    # we don't know what is the order of bond_vec.
    for conf in confs:
      apex_point = conf[0]
      cone_operator = conf[3]
      cone_degree = (conf[2] / 2.0) * numpy.pi/180.0
      end_point_central_vec = None
      end_bond_point = None
      if apex_point == bond[0]:
        apex_point = bond[0]
        end_bond_point = bond[1]
      elif apex_point == bond[1]:
        apex_point = bond[1]
        end_bond_point = bond[0]
      # Now we need to find the second point of the cone vector.
      # It will be an neighbour of the apex_point.
      for nb in apex_point.neighbours:
        if nb == conf[1]:
          end_point_central_vec = nb
          break
      if not end_point_central_vec:
        continue

      central_vec = apex_point.position - end_point_central_vec.position
      bond_vec = end_bond_point.position - apex_point.position
      #print end_bond_point.position, apex_point.position, end_point_central_vec.position
      angle_central_bond = self.angle_between(bond_vec, central_vec)
      #print angle_central_bond, cone_operator, cone_degree
      if cone_operator == '<':
        return angle_central_bond <= cone_degree
      else:
        return angle_central_bond >= cone_degree


  def find_atom_pairs(self, atoms, cutoff):
    """Try to find the atom pairs that should be bonded.

    Args;
      atoms: The list of atoms tuple.
      cutoff: The cutoff distance.

    Returns:
      The list of atoms tuple to be consider as linked.
    """
    cutoff = cutoff ** 2
    valid_atoms = {x.atom_id: x for x in atoms}
    bonds = defaultdict(list)

    l = len(atoms)
    for i in range(l):
      atom_i = atoms[i]
      nbs_of_atom_i = []
      for j in range(i + 1, l, 1):
        atom_j = atoms[j]
        if atom_i.atom_id != atom_j.atom_id:
          atom_pair = self.atom_pairs_settings.get(
              (atom_i.name, atom_i.chain_name, atom_j.name, atom_j.chain_name)
              )
          # Check the intramolecular constraint.
          if not self.local_settings.INTRAMOLECULAR_BONDS:
            if atom_i.chain_idx == atom_j.chain_idx:
              atom_pair = None # Remove the pair, it is in the same
          # So we have atom pair, check the distance constraint and optional
          # check the geometric constraint.
          if atom_pair:
            d_2 = self.calculate_distance(atom_i.position, atom_j.position)
            if d_2 <= cutoff:
              # Geometric constraints.
              valid = True
              for geo_key in self.local_settings.GEOMETRIC_CONSTRAINTS:
                if geo_key == 'cone':
                  ret = self.cone_constraint(atom_pair, (atom_i, atom_j))
                elif geo_key == 'vector':
                  ret = self.vector_constraint(atom_pair, (atom_i, atom_j))
                elif geo_key == 'align':
                  ret = self.align_constraint(atom_pair, (atom_i, atom_j))
                elif geo_key == 'sphere':
                  ret = self.sphere_constraint(atom_pair, (atom_i, atom_j))

                if ret is not None:
                  valid &= ret

                if not valid:  # Well all condition has to valid.
                  break
              if valid:
                nbs_of_atom_i.append((atom_j, d_2))

      # Sort by the distance, so the closest one will be taken.
      nbs_of_atom_i.sort(key=lambda x: x[1])
      for nb in nbs_of_atom_i:
        atom_j = nb[0]
        if ((valid_atoms[atom_i.atom_id].degree <
             self.atoms_def_per_chain[atom_i.chain_name][atom_i.name].max_degree)
            and (valid_atoms[atom_j.atom_id].degree <
                 self.atoms_def_per_chain[atom_j.chain_name][atom_j.name].max_degree
                )):
          atom_pair = self.atom_pairs_settings.get(
              (atom_i.name, atom_i.chain_name, atom_j.name, atom_j.chain_name)
              )
          if atom_pair:
            # Define the new bond.
            bonds[atom_pair].append((atom_i, atom_j))
            atom_i.degree += 1
            atom_j.degree += 1

    return bonds


def find_1_4_pairs(new_bonds, settings):
  """Find the intermolecular pairs.

  Args:
    top_file: The topology file.
    new_bonds: The dictionary with the new bonds, indexed by the bond tuple
    settings: The settings object.

  Returns:
    The dictionary with the keys as the tuple of pair and value with configuration.
  """

  def traverse(atom, s, level=0, prev_atom=None):
    s[level].add(atom)

    level += 1
    if level < 3:
      # Preventing loops.
      nb_list = [x for x in atom.neighbours if not x.cmp_exact(prev_atom)]
      for nb in nb_list:
        traverse(nb, s, level, prev_atom=atom)

  # The 1-4 pairs are the pairs that are separated by 3 atoms.
  #
  new_pairs = {}

  for pair_conf, conf_entry in settings.PAIRS.iteritems():
    bonds = new_bonds[pair_conf]
    for at1, at2 in bonds:
      s1, s2 = defaultdict(set), defaultdict(set)
      # Process atom 1 from the bond. The level 1 because we want to get the information
      # about the molecule id.
      traverse(at1, s1, 0, None)
      traverse(at2, s2, 0, None)

      # The idea is to do the cross product between those two sets but excluding the pairs
      # in which two atoms belongs to the same chain (chain_name and chain_idx are the same)
      # or the pairs that already form the bonds. We based on the neighbours lists of each of
      # the atoms.
      # The cross product is made between the s1[n] x s2[2-n] where n is from <0, 2>. s1 and s2
      # are the dictionaries indexed by the distance from the atom1 and atom2. Those atoms
      # form the new crosslink bond. We assume that the all other pairs are already defined.
      # so only the new pairs are required.
      for l in [0, 1, 2]:
        new_pairs.update({
          (x[0].atom_id, x[1].atom_id): conf_entry
          for x in itertools.product(s1[l], s2[2-l])
          if not (
            (x[0].chain_idx == x[1].chain_idx and x[0].chain_name == x[1].chain_name)
            or ((x[0], x[1]) in bonds or (x[1], x[0]) in bonds)
            )
          })

  return new_pairs

