"""
Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import collections
import networkx as nx


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', help='Config file')
    parser.add_argument('--interactive', action='store_true')
    parser.add_argument('in_file', help='Lammps input file')
    parser.add_argument('out_graph', help='Graph')

    return parser.parse_args()


def lammps_reader(file_name):
    atoms = {}
    bonds = []
    with open(file_name, 'r') as f:
        in_atom_section = False
        in_bond_section = False
        nr_bonds, nr_atoms = 0, 0
        for l in f:
            ll = l.strip()
            if 'atoms' in ll:
                nr_atoms = int(ll.split()[0])
                continue
            elif 'bonds' in ll:
                nr_bonds = int(ll.split()[0])
                continue
            elif 'Atoms' in ll:
                in_atom_section = True
                continue
            elif 'Bonds' in ll:
                in_bond_section = True
                continue
            elif in_atom_section and nr_atoms > 0:
                tmp = ll.split()
                if tmp:
                    atoms[int(tmp[0])] = {'atom_tag': int(tmp[1]), 'atom_type': int(tmp[2])}
                    nr_atoms -= 1
            elif in_bond_section and nr_bonds > 0:
                tmp = ll.split()
                if tmp:
                    bonds.append((int(tmp[1]), int(tmp[2]), int(tmp[3])))
                    nr_bonds -= 1
    return atoms, bonds


def create_graph(atoms, bonds, config):
    atoms_per_mol = collections.defaultdict(list)
    atom_id_mol_id = {}
    g = nx.Graph()
    for at_id, at_val in atoms.iteritems():
        if at_val['atom_type'] in config['ACTIVE_SITES']:
            atoms_per_mol[at_val['atom_tag']].append(at_id)
            atom_id_mol_id[at_id] = at_val['atom_tag']
            g.add_node(at_val['atom_tag'], atom_type=at_val['atom_type'])

    for bond_type, bi, bj  in bonds:
        if bond_type in config['BOND_TYPES']:
            g.add_edge(atom_id_mol_id[bi], atom_id_mol_id[bj], btype=bond_type)

    return g

def main():
    args = _args()
    config = {}
    execfile(args.config, {}, config)

    atoms, bonds = lammps_reader(args.in_file)

    if args.interactive:
        import IPython
        IPython.embed()

if __name__ == '__main__':
    main()
