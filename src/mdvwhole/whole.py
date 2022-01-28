#!/usr/bin/env python
# coding: utf-8

import argparse
import scipy.ndimage
import numpy as np
import mdvoxelsegmentation as mdvseg
import MDAnalysis as mda
import networkx as nx
from time import time
from numba import jit


# Visualization
import open3d as o3d
from pyvis import network as pvnet

# Benchmarking
#%load_ext line_profiler


def make_pcd(atomgroup, color=np.random.random()):
    """
    Returns a pcd with a uniform color as a pcd.
    """
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(atomgroup.positions)
    pcd.paint_uniform_color(color)
    return pcd


def draw_atomgroup(atomgroup, color=np.random.random(3)):
    """
    Draws an atomgroup as a pointcloud using o3d.
    """
    pcd = make_pcd(atomgroup, color)
    vis = o3d.visualization.draw_geometries([pcd])


def draw_grid(grid):
    """
    Draws a boolean or integer grid using o3d. Can also be used
    to visualize the labels in voxel format.
    """
    labels = np.unique(grid)
    
    # Create some random colors and prevent complete white.
    colors = np.random.random((len(labels), 3))
    colors[colors >= 0.95] = 0.95
    
    pcds = []
    for idx, label in enumerate(labels):
        if label == 0:
            continue
        points = np.array(np.where(grid == label)).T
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(points)
        pcd.paint_uniform_color(colors[idx])
        pcds.append(pcd)
    o3d.visualization.draw_geometries(pcds)


@jit(nopython=True)
def find_neighbor_indexes(target_indexes, labels, neighbor_mask):
    """
    Returns the indexes of the 27 neighbors and the dim shift in two
    seperate arrays of (3,3,3) int32.
    """
    neighbor_indexes = np.zeros((3,3,3,3), dtype='int32')
    neighbor_indexes += np.array(target_indexes, dtype='int32')
    neighbor_indexes += neighbor_mask
    neighbor_dim_shifts, neighbor_indexes = np.divmod(neighbor_indexes, np.array(labels.shape))
    return neighbor_indexes, neighbor_dim_shifts


@jit(nopython=True)
def find_bridges(target_indexes, labels, neighbor_mask):
    """
    Returns the bridges with the center label, its connecting labels 
    and the corresponding dimension shifts.
    """
    target_label = labels[target_indexes]
    bridges = np.zeros((27), dtype='int32')
    bridge_shifts = np.zeros((27,3), dtype='int32')
    if target_label != 0:
        indexes, dim_shifts = find_neighbor_indexes(target_indexes, labels, neighbor_mask)
        counter = 0
        for idxx, x in enumerate(indexes):
            for idxy, y in enumerate(x):
                for idxz, z in enumerate(y):
                    bridges[counter] = labels[z[0], z[1], z[2]]
                    bridge_shifts[counter] = dim_shifts[idxx, idxy, idxz]
                    counter += 1
    return target_label, bridges, bridge_shifts


@jit(nopython=True)
def find_all_bridges(labels, neighbor_mask):
    """
    Returns all bridges for all voxels. 
    
    Returning a list of the voxel bridges with the center label, 
    its connecting labels and the corresponding dimension shifts.
    """
    voxel_amount = labels.shape[0]*labels.shape[1]*labels.shape[2]
    all_target_labels = np.zeros((voxel_amount), dtype='int32')
    all_pairs = np.zeros((voxel_amount, 27))
    all_pair_shifts = np.zeros((voxel_amount, 27, 3))
    counter = 0
    for idxx, x in enumerate(labels):
        for idxy, y in enumerate(x):
            for idxz, z in enumerate(y):
                all_target_labels[counter], all_pairs[counter], all_pair_shifts[counter] = find_bridges(
                    (idxx, idxy, idxz), labels, neighbor_mask)
                counter += 1
    return all_target_labels, all_pairs, all_pair_shifts


def create_edge_voxel_array(labels):
    """
    Returns all edge voxels (this could be cheaper but I do not care for now).
    """
    edge_mask = np.zeros(labels.shape, dtype='int32')
    edge_mask[-1, :, :] = 1
    edge_mask[:, -1, :] = 1
    edge_mask[:, :, -1] = 1
    return np.array(np.where(edge_mask == 1), dtype='int32').T


@jit(nopython=True)
def find_all_bridges_from_array(query_array, labels, neighbor_mask):
    """
    Returns all bridges for all edge voxels (one sided).
    
    Returning a list of the voxel pairs with the center label, 
    its connecting labels and the corresponding dimensifts.
    """
    voxel_amount = len(query_array)
    all_target_labels = np.zeros((voxel_amount), dtype='int32')
    all_bridges = np.zeros((voxel_amount, 27))
    all_bridge_shifts = np.zeros((voxel_amount, 27, 3))
    for idx, query in enumerate(query_array):
        all_target_labels[idx], all_bridges[idx], all_bridge_shifts[idx] = find_bridges(
                    (query[0], query[1], query[2]), labels, neighbor_mask)
    return all_target_labels, all_bridges, all_bridge_shifts


@jit(nopython=True)
def create_graph_input_single(label, bridges, dim_shifts):
    """
    Creates an array to generate the graph from the connectivity for one voxel.
    """
    bridges = bridges.astype('int32')
    dim_shifts = dim_shifts.astype('int32')
    if label != 0 and len(np.unique(bridges)) > 2:
        non_self = np.where((bridges != 0) & (bridges != label))[0].astype('int32')
        out = np.zeros((len(non_self), 5), dtype='int32')
        for idx, single_pointer in enumerate(non_self):
            out[idx] = (label, bridges[single_pointer], dim_shifts[single_pointer][0], 
            dim_shifts[single_pointer][1], dim_shifts[single_pointer][2])
        return out
    else:
        return None


def filter_largest_bridges(graph_input):
    """
    Creates a sorted array to generate the graph from the connectivity for all voxels.
    
    the bridges are sorted from high to low occurance.
    
    !!!TODO!!! (or check)
    I think it would be better if I would only return the most bridged dimension.
    This is eventually what happens later on in the graph path finding. This is not
    clear in the code at all, but for some reason the first is picked in the multigraph,
    maybe this actually breaks on the (double) diagonal, for here the node count would
    actually be shorter. Thus, remove all but the largest connection migh be more
    predictable, on the other hand, would we miss something is we only kept the largest?
    
    label, target_label, x, y, z
    """
    sorted_graph_input = np.copy(graph_input)
    # Use the minus sign to invert the search to descending.
    sorted_graph_input = sorted_graph_input[np.argsort(-sorted_graph_input[:, 5])]
    return sorted_graph_input


def create_graph_input_all(labels, bridges, dim_shifts):
    """
    Creates an array to generate the graph from the connectivity for all voxels.
    
    label, target_label, x, y, z
    """
    all_arrays = []
    for idx, _ in enumerate(labels):
        temp_out = create_graph_input_single(labels[idx], bridges[idx], dim_shifts[idx])
        if temp_out is not None:
            all_arrays.append(temp_out)
    try:
        all_arrays = np.vstack(all_arrays)
    except ValueError:
        return np.zeros((2,26))
    # We need to check for the largest PBC connection per label.
    temp_out = np.unique(all_arrays, axis=0, return_counts=True)
    out = np.zeros((temp_out[1].shape[0], 6))
    out[:,:5] = temp_out[0]
    out[:, 5] = temp_out[1]
    out = filter_largest_bridges(out)
    return out


class Voxels():
    def __init__(self, atomgroup, resolution=1, hyperres=False):
        self.atomgroup = atomgroup
        # TODO set the atoms back in their PBC using a method
        #  which works for all triclinic boxes.
        #pbc = mdv.clustering. self.atomgroup.dimensions
        self.atomgroup.positions %= self.atomgroup.dimensions[:3]
        self._voxelate(resolution, hyperres)
        self._neighbor_mask = self._generate_neighbor_mask()
        
    def _generate_neighbor_mask(self):
        """
        Returns the offset for a 27 neighbor connection.
        
        #TODO Update this mask and all functions using it
        so that it also allows for a 1/8 top right view.
        This due to the fact that every contact is symmetric
        in every dimension. This makes the whol a lot faster.
        However, the select neighbors is nice anyway to have
        so I think I will keep it like this and add the 
        other option for speed ups in certain scenarios
        with symmetries (like this algorithm).
        """
        neighbor_mask = np.zeros((3,3,3,3), dtype='int32')
        neighbor_mask[1,0,0] = ( 0, -1, -1)
        neighbor_mask[2,0,0] = ( 1, -1, -1)

        neighbor_mask[0,1,0] = (-1,  0, -1)
        neighbor_mask[0,2,0] = (-1,  1, -1)

        neighbor_mask[0,0,1] = (-1, -1,  0)
        neighbor_mask[0,0,2] = (-1, -1,  1)

        neighbor_mask[1,1,0] = ( 0,  0, -1)
        neighbor_mask[0,1,1] = (-1,  0,  0)
        neighbor_mask[1,0,1] = ( 0, -1,  0)

        neighbor_mask[2,2,0] = ( 1,  1, -1)
        neighbor_mask[0,2,2] = (-1,  1,  1)
        neighbor_mask[2,0,2] = ( 1, -1,  1)

        neighbor_mask[1,1,2] = ( 0,  0,  1)
        neighbor_mask[2,1,1] = ( 1,  0,  0)
        neighbor_mask[1,2,1] = ( 0,  1,  0)

        neighbor_mask[2,2,1] = ( 1,  1,  0)
        neighbor_mask[1,2,2] = ( 0,  1,  1)
        neighbor_mask[2,1,2] = ( 1,  0,  1)

        neighbor_mask[0,1,2] = (-1,  0,  1)
        neighbor_mask[1,2,0] = ( 0,  1, -1)
        neighbor_mask[2,1,0] = ( 1,  0, -1)

        neighbor_mask[1,0,2] = ( 0, -1,  1)
        neighbor_mask[0,2,1] = (-1,  1,  0)
        neighbor_mask[2,0,1] = ( 1, -1,  0)

        neighbor_mask[0,0,0] = (-1, -1, -1)
        neighbor_mask[1,1,1] = ( 0,  0,  0)
        neighbor_mask[2,2,2] = ( 1,  1,  1)
        return neighbor_mask
        
    def _voxelate(self, resolution=1, hyperres=False):
        """
        Sets the boolean grid (self.grid) using the atomgroup. This
        respects all triclinic box types. It also sets the used hyperres
        and resolution values for later reference.
        """
        self.resolution = resolution
        self.hyperres = hyperres
        self.grid, self.voxel2atom, self.nbox = mdvseg.clustering.gen_explicit_matrix(
            self.atomgroup, resolution, hyperres)
        
    def get_voxel(self, coordinates=(0,0,0)):
        """
        Returns the atomgroup of all atoms in the specified voxel (xyz).
        """
        atomids = self.voxel2atom[coordinates]
        atomgroup = self.atomgroup[atomids]
        return atomgroup
    
    def get_voxels(self, voxels=((0,0,0), (0,0,1))):
        """
        Returns the atomgroup of all the atoms in the specified list 
        of voxels (xyz).
        """
        atomgroup = mdvseg.clustering.voxels2atomgroup(
            voxels, self.voxel2atom, self.atomgroup)
        return atomgroup
    
    def find_neighbor_indexes(self, target_indexes):
        """
        Returns the indexes of the neighboring voxels.
        """
        return find_neighbor_indexes(target_indexes, self.labels, self._neighbor_mask)
    
    def _label_count(self):
        """
        Returns a dictionary with the labels as keys and the amount of 
        voxels in the label as count.
        """
        unique, counts = np.unique(self.labels, return_counts=True)
        return dict(zip(unique, counts))
    
    def label(self, neighbor_mask=np.ones((3,3,3))):
        """
        Sets (self.labels) and returns the label masked grid. The
        neighbor_mask specifies what is considered a contact.
        """
        self.neighbor_mask = neighbor_mask
        labels = scipy.ndimage.label(self.grid, neighbor_mask)
        self.labels = labels[0]
        self.label_count = self._label_count()
        return self.labels
    
    def __checkNset_labels(self):
        """
        Checks if the labels have been generated, if not, generate them using the default
        neighbormask of 26 neighbors.
        """
        if not hasattr(self, 'labels'):
            print('The default 26 neighbormask has been used for labeling.\n'
                  'Running Voxel.label() allows for specification of a custom\n'
                  'neighbormask.\n')
            self.label()
    
    def get_label_voxels(self, labelID):
        """
        Returns the list of voxels (xyz) which make up the specified label.
        """
        # If not self.labels, create them.
        self.__checkNset_labels()
            
        # Get the voxels.
        voxels = np.array(np.where(self.labels == labelID)).T
        voxels = [tuple(voxelID) for voxelID in voxels]
        return voxels
    
    def get_label(self, labelID):
        """
        Returns the atomgroup of all atoms which lie in a voxel with the
        specified labelID (int).
        """
        # If not self.labels, create them.
        self.__checkNset_labels()
        
        # Get the atomgroup.
        voxels = self.get_label_voxels(labelID)
        atomgroup = self.get_voxels(voxels)
        return atomgroup
    
    def get_labels(self):
        """
        Sets (self.label_atomgroups) and returns a dictionary with the 
        labelIDs as the keys and the label atomgroups as values.
        """
        # If not self.labels, create them.
        self.__checkNset_labels()
        
        # Get the labels atomgroups. 
        atomgroups = []
        labelIDs = list(self.label_count.keys())
        for label in labelIDs:
            atomgroups.append(self.get_label(label))
        self.label_atomgroups = dict(zip(labelIDs, atomgroups))
        return self.label_atomgroups
    
    def draw_label(self, label, color=np.random.random(3)):
        """
        Draws all the atoms in the label using o3d. A color can be
        specified.
        """
        # If not labels, create them.
        self.__checkNset_labels()
        
        # Draw the label.
        label_atomgroup = self.get_label(label)
        print(f'Drawing label {label} with {len(label_atomgroup)} atoms.')
        draw_atomgroup(label_atomgroup, color)
        
    def draw_labels(self):
        """
        Draws all labels with a random color using o3d.
        
        #TODO Use a seaborn cycle for colors with specified spread.
        """
        # If not self.label_atomgroups create them.
        if not hasattr(self, 'label_atomgroups'):
            self.get_labels()
        
        # Perform the drawing logic.
        pcds = []
        # Create some random colors and prevent complete white.
        colors = np.random.random((len(self.label_count), 3))
        colors[colors >= 0.95] = 0.95
        # Create the pointclouds
        for idx, atomgroup in enumerate(self.label_atomgroups.values()):
            pcds.append(make_pcd(atomgroup, colors[idx]))
        # Draw the pointclouds.
        o3d.visualization.draw_geometries(pcds)
        
    def draw_grid(self):
        """
        Draws the grid using o3d.
        """
        draw_grid(self.grid)
    
    def draw_labels_grid(self):
        """
        Draws the labels using o3d.
        """
        # If not labels, create them.
        self.__checkNset_labels()   
        draw_grid(self.labels)
    

class Bridges():
    def __init__(self, labels, neighbor_mask):
        self.labels = labels
        self.labelIDs = np.unique(labels)
        self._neighbor_mask = neighbor_mask
        self.bridges = self._find_all_bridges()
        filter_largest_bridges(self.bridges)
        self._make_graph()
        self._find_subgraphs()
        
    def _find_bridges(self, target_indexes):
        """
        Returns an array of all bridges for the target voxel.
        
        Bridges are annotated with their associated dimension
        shift.
        
        np.array():
            label, connecting_label, x, y, z
        """
        return find_bridges(target_indexes, self.labels, self._neighbor_mask)
    
    def _find_all_bridges(self):
        """
        Returns an array of all bridges.
        
        Bridges are annotated with their associated dimension
        shift.
        
        np.array():
            label, connecting_label, x, y, z
        """
        indexes = create_edge_voxel_array(self.labels)
        #indexes = np.array(np.where(self.labels[indexes.T] >= 1), dtype='int32').T
        results = find_all_bridges_from_array(indexes, self.labels, self._neighbor_mask)
        return create_graph_input_all(*results) 
    
    def _make_graph(self):
        """
        Sets and returns a directed multigraph of the bridge connectivity. 
        The dimension, direction and connectivity are stored in the edge 
        as an special attribute 'dimension' and the connectivity is the
        weight of the edge. 
        """
        G = nx.MultiDiGraph()
        G.add_nodes_from(range(len(self.labelIDs)))
        for bridge in self.bridges:
            G.add_edge(int(bridge[0]), int(bridge[1]), label=f'{ bridge[2:5]}', value=int(bridge[5]), cost= tuple(bridge[2:5]))
            G.add_edge(int(bridge[1]), int(bridge[0]), label=f'{-bridge[2:5]}', value=int(bridge[5]), cost= tuple(-bridge[2:5]))
        self.graph = G
        return G
    
    def _find_subgraphs(self):
        """
        Sets and returns the subgraphs in the bridges as list of set(labelIDs)
        for every subgraph in the bridges.
        """
        # make an undirected copy of the digraph
        UG = self.graph.to_undirected()
    
        # extract subgraphs
        sub_graphs = list(nx.connected_components(UG))
        graph_ids = range(len(sub_graphs))
        self.sub_graphs_nodes = dict(zip(graph_ids, sub_graphs))
        return self.sub_graphs_nodes
        
    
    def draw_graph(self, name='out.html', height='600px', width='600px'):
        """
        Draws the bridges graph using pvnet. The nodes represent the labelIDs.
        The width of the edges represents the amount of bridges between two 
        segements. The edge label is the dimension of the connection.
        """
        g = self.graph.copy() # some attributes added to nodes
        net = pvnet.Network(notebook=True, directed=True, height=height, width=width)
        opts = '''
            var options = {
              "physics": {
                "forceAtlas2Based": {
                  "gravitationalConstant": -200,
                  "centralGravity": 0.11,
                  "springLength": 100,
                  "springConstant": 0.09,
                  "avoidOverlap": 1
                },
                "minVelocity": 0.75,
                "solver": "forceAtlas2Based",
                "timestep": 0.22
              }
            }
        '''

        net.set_options(opts)
        # uncomment this to play with layout
        #net.show_buttons(filter_=['physics'])
        net.from_nx(g)
        return net.show(name)


class Whole():
    def __init__(self, atomgroup, resolution=1, hyperres=False, 
                 neighbor_mask=np.ones((3,3,3))):
        self.voxels = Voxels(atomgroup, resolution, hyperres)
        self.voxels.label(neighbor_mask=neighbor_mask)
        self.bridges = Bridges(self.voxels.labels, self.voxels._neighbor_mask)
        self._find_biggest_labels()
        self._find_all_paths()
        self._shift_all_atoms()
        
    def _find_biggest_labels(self):
        """
        Sets (self.biggest_labels) and returns a dictionary with the 
        graphIDs as the keys and the biggest label in the subgraph 
        as the value. 
        """
        biggest_labels = dict()
        for graphID in self.bridges.sub_graphs_nodes.keys():
            graph_sizes= dict()
            for label in self.bridges.sub_graphs_nodes[graphID]:
                if label == 0:
                    continue
                label_size = len(self.voxels.get_label_voxels(label))
                graph_sizes[label] = label_size
            try:
                biggest_label = max(graph_sizes, key=graph_sizes.get)
                biggest_labels[graphID] = biggest_label
            except ValueError:
                pass
        self.biggest_labels = biggest_labels
        return biggest_labels
    
    def _find_paths(self, sub_graph_id, center_label=False):
        """
        Rerturns all shortest paths to the central label for a sub_graph.
        """
        sub_graph = self.bridges.graph.subgraph(self.bridges.sub_graphs_nodes[sub_graph_id])
        if center_label is False:
            center_label = self.biggest_labels[sub_graph_id]
        shortest_paths = nx.algorithms.shortest_paths.generic.shortest_path(sub_graph, center_label)
        return shortest_paths
    
    def _find_all_paths(self):
        """
        Sets (self.all_paths) and returns all shortest paths for 
        all subgraphs as dict_sub_graph(dict_paths()).
        """
        all_paths = dict()
        for sub_graph_id in self.bridges.sub_graphs_nodes.keys():
            # This is needed for the 0 label. It is a bit strange and 
            #  does not have an entry in the biggest_labels.
            try:
                all_paths[sub_graph_id] = self._find_paths(sub_graph_id)
            except KeyError:
                continue
        self.all_paths = all_paths
        return all_paths

    def _calculate_shift(self, sub_graph_id, label):
        """
        Returns the shift in box dimensions for the label in the sub_graph.
        """
        path = self.all_paths[sub_graph_id][label]
        counter = 0
        path_sum = np.zeros(3)
        while counter < len(path)-1:
                # Take the first edge in the dict. MultiGraph allows for multiple 
                #  paths which could happen with multiple dimensions. However,
                #  we, for now, always select the first. This is rather tricky
                #  to wrap my head around and maybe this requires extra thought.
                #  For now lets hope it all works fine :D
                edge_data = self.bridges.graph.get_edge_data(path[counter], path[counter+1])[0]
                path_sum += edge_data['cost']
                counter += 1
        return path_sum
    
    def _shift_atoms(self, label, shift):
        """
        Returns the atomgroup with the label shifted.
        """
        shifted_atomgroup = self.voxels.atomgroup
        label_atomgroup = self.voxels.get_label(label)
        box_dimensions = shifted_atomgroup.dimensions
        shifted_atomgroup.intersection(label_atomgroup).positions += box_dimensions[:3] * shift

        return shifted_atomgroup
    
    def _shift_all_atoms(self):
        """
        Return the shifted atomgroup. Shifts all the atoms in all sub_graphs.

        #TODO find out how to make this all happen in a copy. This prevents the 
        #  multiple shifting of positions if the method is ran more than once.
        #  also it is neat to leave the original data unaffected? This would
        #  increase the amount of memory used.
        """
        for sub_graph_id in self.bridges.sub_graphs_nodes.keys():
            for label in self.bridges.sub_graphs_nodes[sub_graph_id]:
                if label == 0:
                    continue
                shift = self._calculate_shift(sub_graph_id, label)
                self._shift_atoms(label, shift)
        return self.voxels.atomgroup
    
    def write_atomgroup(self, name='whole.gro'):
        """
        Writes the whole atomgroup.
        """
        self.voxels.atomgroup.write(name)


def test(gro, xtc, frame, selection='not resname W WF ION', resolution=1, name='whole.gro'):
    # Load the universe.
    u = mda.Universe(gro, xtc)
    atomgroup = u.select_atoms(selection)
    #u.trajectory[1200]
    #u.trajectory[200]
    u.trajectory[frame]
    
    draw_atomgroup(atomgroup)

    # Create the voxel instance for testing.
    voxels = Voxels(atomgroup, resolution=resolution)
    voxels.label(neighbor_mask=np.ones((3,3,3)))
    voxels.draw_labels()

    # Create the whole instance for testing.
    start_time = time()
    whole = Whole(atomgroup, resolution=resolution)
    total_time = time() - start_time
    print(f'It took {total_time:.4f} seconds to make whole.')
    graph = whole.bridges.draw_graph()
    whole.voxels.draw_labels()
    
    # Draw the final atomgroup without and with labels
    #  in the make_whole state
    start_time = time()
    whole.write_atomgroup(name)
    total_time = time() - start_time
    print(f'It took {total_time:.4f} seconds to write whole.')
    draw_atomgroup(whole.voxels.atomgroup)
    #TODO make it such so I can draw the labeling in the new make_whole.
    #  This is not working for the labeling is a voxel thing, the voxels 
    #  are filled with PBC!!! We should be able to turn this off...
    
    return graph


class MDAWhole():
    """
    A MDAnalysis compatible Whole class which can be used
    as a just in time transformation or to make the complete
    trajectory whole and write the output.
    """
    def __init__(self, atomgroup, resolution=1):
        """
        Sets the atomgroup.
        """
        self.atomgroup = atomgroup
        self.resolution = resolution
    def __call__(self, ts):
        """
        Used in the MDA transformations logic for making a frame
        whole just in time.
        """
        Whole(self.atomgroup, resolution=self.resolution)
        return ts
    @classmethod
    def whole_traj(cls, atomgroup, resolution=1, out='test.xtc'):
        """
        Makes every frame whole, writes the xtc and returns the 
        whole atomgroup.
        """
        u = atomgroup.universe
        total_frames = len(u.trajectory)
        start_time = time()
        u.trajectory.add_transformations(MDAWhole(atomgroup, resolution))
        with mda.Writer(out, atomgroup.n_atoms) as W:
            for frame in u.trajectory:
                # Making a good estimate of the remainig time
                current_frame_id = u.trajectory.frame
                frame_time = time()
                projected_total_time = ((total_frames/(current_frame_id+0.0001))) * (frame_time-start_time)
                time_left = (1 - current_frame_id/total_frames) * projected_total_time
                # Some time left printing logic (seconds, minutes, hours)
                if time_left <= 60:
                    print(f'\rFrame {current_frame_id}/{total_frames} {time_left:.2f} seconds remaining.     ', end='')
                elif time_left < 3600:
                    time_left /= 60
                    print(f'\rFrame {current_frame_id}/{total_frames} {time_left:.2f} minutes remaining.     ', end='')
                else:
                    time_left /= 3600
                    print(f'\rFrame {current_frame_id}/{total_frames} {time_left:.2f} hours remaining.     ', end='')
                # Doing the actual calculation
                W.write(atomgroup)
        print(f'\rDone, the whole thing took {(time()-start_time)/60:.2f} minutes.              ')
        return atomgroup


def test(gro, xtc, frame, selection='not resname W WF ION', resolution=1, name='whole.gro'):
    # Load the universe.
    u = mda.Universe(gro, xtc)
    atomgroup = u.select_atoms(selection)
    #u.trajectory[1200]
    #u.trajectory[200]
    u.trajectory[frame]
    
    draw_atomgroup(atomgroup)

    # Create the voxel instance for testing.
    voxels = Voxels(atomgroup, resolution=resolution)
    voxels.label(neighbor_mask=np.ones((3,3,3)))
    voxels.draw_labels()

    # Create the whole instance for testing.
    start_time = time()
    whole = Whole(atomgroup, resolution=resolution)
    total_time = time() - start_time
    print(f'It took {total_time:.4f} seconds to make whole.')
    graph = whole.bridges.draw_graph()
    whole.voxels.draw_labels()
    
    # Draw the final atomgroup without and with labels
    #  in the make_whole state
    start_time = time()
    whole.write_atomgroup(name)
    total_time = time() - start_time
    print(f'It took {total_time:.4f} seconds to write whole.')
    draw_atomgroup(whole.voxels.atomgroup)
    #TODO make it such so I can draw the labeling in the new make_whole.
    #  This is not working for the labeling is a voxel thing, the voxels 
    #  are filled with PBC!!! We should be able to turn this off...
    
    return graph


def whole_traj(atomgroup, resolution=1, out='whole.xtc'):
    """
    Makes every frame whole and write the xtc and returns the whole atomgroup.
    """
    u = atomgroup.universe
    total_frames = len(u.trajectory)
    start_time = time()
    with mda.Writer(out, atomgroup.n_atoms) as W:
        for frame in u.trajectory:
            # Making a good estimate of the remainig time
            current_frame_id = atomgroup.universe.trajectory.frame
            frame_time = time()
            projected_total_time = ((total_frames/(current_frame_id+0.0001))) * (frame_time-start_time)
            time_left = (1 - current_frame_id/total_frames) * projected_total_time
            # Some time left printing logic (seconds, minutes, hours)
            if time_left <= 60:
                print(f'\rFrame {current_frame_id}/{total_frames} {time_left:.2f} seconds remaining.     ', end='')
            elif time_left < 3600:
                time_left /= 60
                print(f'\rFrame {current_frame_id}/{total_frames} {time_left:.2f} minutes remaining.     ', end='')
            else:
                time_left /= 3600
                print(f'\rFrame {current_frame_id}/{total_frames} {time_left:.2f} hours remaining.     ', end='')
            # Doing the actual calculation
            whole = Whole(atomgroup, resolution=resolution)
            W.write(whole.voxels.atomgroup)     
    print(f'\rDone, the whole thing took {(time()-start_time)/60:.2f} minutes.              ')
    return atomgroup


def read_arguments():
    """
    Parses the input arguments from the command line.
    Returns
    -------
    args = NameSpace
    """
    # Generating the argparser object.
    parser = argparse.ArgumentParser(add_help=False)
    required_grp = parser.add_argument_group(title='required arguments')
    optional_grp = parser.add_argument_group(title='optional arguments')
    # REQUIRED
    required_grp.add_argument(
        '-f', '--reference', nargs='?', required=True, type=str,
        help='an MDAnalysis compatible coordinate file including atom names, resids and positions (e.g. GRO, TPR or PDB)',
        )
    required_grp.add_argument(
        '-x', '--trajectory', nargs='?', required=True, type=str,
        help='an MDAnalysis compatible trajectory file (e.g. XTC)',
        )
    
    # OPIONAL
    optional_grp.add_argument(
        '-sel', '--selection', nargs='?', default='not resname W WF ION', 
        type = str,
        help='the selection excluding the solvent',
        )    
    optional_grp.add_argument(
        '-res', '--resolution', nargs='?', default=1, type=float,
        help='the binning resolution in the same units as the reference file (default=1 nm)',
        )
    optional_grp.add_argument(
        '-o', '--out_file', nargs='?', default='whole.xtc', type=str,
        help='the path for writing the whole (e.g. XTC, GRO) (default=whole.xtc)',
        )
    optional_grp.add_argument(
    '-h', '--help', action="help",
    help='show this help message and exit',
    )
    # parse the arguments into the name space
    args = parser.parse_args()
    return args


def main():
    """
    Makes the trajectory whole and writes it.
    """
    args = read_arguments()
    print(args)
    # Loading trajectory.
    u = mda.Universe(args.reference, args.trajectory, in_memory=False)
    atomgroup = u.select_atoms(args.selection)
    # Make the complete trajectory whole and write it.
    MDAWhole.whole_traj(atomgroup, out=args.out_file, resolution=args.resolution)


if __name__ == "__main__":
    main()


# # Performance
# The main performance hit is the creation of the voxels2atoms dict. I could spent time on only checking half the neighbors due to symmetry, but it doesn't seem to be worth the effort. I think we could try no using a dict but the array itself, using np.unique and np.where(unique_value) instead. I am not sure if this will be faster in the end as you have to perform the np.where(unique_value) check for every unique value. This might turn out to be more expensive as the array gets bigger?
# 
# # To be added
# * A decent argparser.
# * Shifting non-selected particles with the shift of the label they are closest to (associative whole).
# * Submasking molecules e.g. segmenting only the tails, but projecting the shift on the whole lipid (sub whole). Sub whole could also work using the associate whole method, but I think it better is to make a separate function. Since it is cheaper if you do not have to do the distance check to know which label you belong to, as you are the same molecule!
