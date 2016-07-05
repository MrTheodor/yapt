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

Network analysis - basic graphs analysis."""

import networkx


def prepare_graph(input_file, topology_file):
  """Prepares graph object from the input file and the topology file.
  
  Args:
    input_file: The input File object.
    topology_file: The topology object.
 
 Returns:
  Return the graph.
  """

  if not input_file.opened:
    raise Exception('The input file has to be read first!')
  
  if not topology_file.opened:
    raise Exception('The topology file has to be read first!')

  G = networkx.Graph()
  
  # Creates the nodes.
  for node_id, node_data in input_file.data.iteritems():
    G.add_node(node_id, node_data.get_dict())
  
  # Creates edges.
  for bond in topology_file.data['bonds']:
    G.add_edge(*bond)

  return G
