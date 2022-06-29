#!/usr/bin/env python
# coding: utf-8

import argparse
import scipy.ndimage
import numpy as np
import mdvoxelsegmentation as mdvseg
import MDAnalysis as mda
import networkx as nx
from time import time
from MDAnalysis import transformations


# Visualization
#import open3d as o3d
#from pyvis import network as pvnet

# Benchmarking
#%load_ext line_profiler


def dim2lattice(x, y, z, alpha=90, beta=90, gamma=90):
    """Convert dimensions (lengths/angles) to lattice matrix"""
    cosa = np.cos( np.pi * alpha / 180 )
    cosb = np.cos( np.pi * beta / 180 )
    cosg = np.cos( np.pi * gamma / 180 )
    sing = np.sin( np.pi * gamma / 180 )

    zx = z * cosb
    zy = z * ( cosa - cosb * cosg ) / sing
    zz = np.sqrt( z**2 - zx**2 - zy**2 )

    return np.array([x, 0, 0, y * cosg, y * sing, 0, zx, zy, zz]).reshape((3,3))


# =============================================================================
# def make_pcd(atomgroup, color=np.random.random()):
#     """
#     Returns a pcd with a uniform color as a pcd.
#     """
#     pcd = o3d.geometry.PointCloud()
#     pcd.points = o3d.utility.Vector3dVector(atomgroup.positions)
#     pcd.paint_uniform_color(color)
#     return pcd
# 
# 
# def draw_atomgroup(atomgroup, color=np.random.random(3)):
#     """
#     Draws an atomgroup as a pointcloud using o3d.
#     """
#     pcd = make_pcd(atomgroup, color)
#     o3d.visualization.draw_geometries([pcd])
# 
# 
# def draw_grid(grid):
#     """
#     Draws a boolean or integer grid using o3d. Can also be used
#     to visualize the labels in voxel format.
#     """
#     labels = np.unique(grid)
#     
#     # Create some random colors and prevent complete white.
#     colors = np.random.random((len(labels), 3))
#     colors[colors >= 0.95] = 0.95
#     
#     pcds = []
#     for idx, label in enumerate(labels):
#         if label == 0:
#             continue
#         points = np.array(np.where(grid == label)).T
#         pcd = o3d.geometry.PointCloud()
#         pcd.points = o3d.utility.Vector3dVector(points)
#         pcd.paint_uniform_color(colors[idx])
#         pcds.append(pcd)
#     o3d.visualization.draw_geometries(pcds)
# 
# =============================================================================

class Voxels():
    def __init__(self, atomgroup, resolution=1, hyperres=False):
        # TODO set the atoms back in their PBC using a method
        #  which works for all triclinic boxes.
        self.atomgroup = atomgroup
        self._voxelate(resolution, hyperres)
        # self.neighbor_mask = np.mgrid[-1:2, -1:2, -1:2].T
        
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
    
# =============================================================================
#     def draw_label(self, label, color=np.random.random(3)):
#         """
#         Draws all the atoms in the label using o3d. A color can be
#         specified.
#         """
#         # If not labels, create them.
#         self.__checkNset_labels()
#         
#         # Draw the label.
#         label_atomgroup = self.get_label(label)
#         print(f'Drawing label {label} with {len(label_atomgroup)} atoms.')
#         draw_atomgroup(label_atomgroup, color)
# =============================================================================
        
# =============================================================================
#     def draw_labels(self):
#         """
#         Draws all labels with a random color using o3d.
#         
#         #TODO Use a seaborn cycle for colors with specified spread.
#         """
#         # If not self.label_atomgroups create them.
#         if not hasattr(self, 'label_atomgroups'):
#             self.get_labels()
#         
#         # Perform the drawing logic.
#         pcds = []
#         # Create some random colors and prevent complete white.
#         colors = np.random.random((len(self.label_count), 3))
#         colors[colors >= 0.95] = 0.95
#         # Create the pointclouds
#         for idx, atomgroup in enumerate(self.label_atomgroups.values()):
#             pcds.append(make_pcd(atomgroup, colors[idx]))
#         # Draw the pointclouds.
#         o3d.visualization.draw_geometries(pcds)
# =============================================================================
        
# =============================================================================
#     def draw_grid(self):
#         """
#         Draws the grid using o3d.
#         """
#         draw_grid(self.grid)
# =============================================================================
    
# =============================================================================
#     def draw_labels_grid(self):
#         """
#         Draws the labels using o3d.
#         """
#         # If not labels, create them.
#         self.__checkNset_labels()   
#         draw_grid(self.labels)
# =============================================================================
    

class Bridges:
    def __init__(self, labelarray, vbox):
        # Set the mask for the faces of the array (given voxel distance d)
        d = 1
        vx, vy, vz = labelarray.shape
        grid = np.mgrid[:vx, :vy, :vz].transpose((1,2,3,0))
        faces = ~np.pad(np.ones((vx-2*d, vy-2*d, vz-2*d) ,dtype=bool), d)

        # Get the labels and voxelindices for non-empty face voxels
        facelabels = labelarray[faces]
        labeled = facelabels > 0
        facevoxels = grid[faces][labeled]
        facelabels = facelabels[labeled]

        # Group the voxels per label
        order = facelabels.argsort()
        facevoxels = facevoxels[order]
        facelabels = facelabels[order]
        labels, counts = np.unique(facelabels, return_counts=True)
        perlabel = np.split(facevoxels, counts[:-1].cumsum())

        # Restrict the shifts 
        up = [] # Don't shift these
        lo = [] # Shift only these
        for group in perlabel:
            up.append(group[np.any(group >= np.diagonal(vbox) - d, axis=1)])
            lo.append(group[np.any(group < d, axis=1)])
        
        # These are the 7 positive lattice shifts
        shifts = np.mgrid[:2,:2,:2].T.reshape((-1,3))[1:]

        # For every combination and every shift determine the bridges
        bridges = np.zeros((len(labels), len(labels), len(shifts)), dtype=int)
        for i, A in enumerate(lo):
            for s, shift in enumerate(shifts):
                shifted = A.copy()
                for k in range(3):
                    if shift[k]:
                        shifted = shifted[shifted[:, k] <= d]
                if len(shifted) == 0:
                    continue
                shifted = shifted + (shift @ vbox)
                for j, B in enumerate(up):
                    bridges[i, j, s] = np.sum(np.abs(shifted[:, None] - B[None, :]).max(axis=2) <= 1)

        # make_graph
        G = nx.MultiDiGraph()
        G.add_nodes_from(labels)

        for a, b, s in zip(*np.where(bridges)):
            shift = shifts[s]
            G.add_edge(labels[a], labels[b], label=str(-shift), value=bridges[a,b,s], cost=-shift)
            G.add_edge(labels[b], labels[a], label=str(shift), value=bridges[a,b,s], cost=shift)

        self.graph = G
        self._find_subgraphs()
        
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
    
# =============================================================================
#     def draw_graph(self, name='out.html', height='600px', width='600px'):
#         """
#         Draws the bridges graph using pvnet. The nodes represent the labelIDs.
#         The width of the edges represents the amount of bridges between two 
#         segements. The edge label is the dimension of the connection.
#         """
#         g = self.graph.copy() # some attributes added to nodes
#         net = pvnet.Network(notebook=True, directed=True, height=height, width=width)
#         opts = '''
#             var options = {
#               "physics": {
#                 "forceAtlas2Based": {
#                   "gravitationalConstant": -200,
#                   "centralGravity": 0.11,
#                   "springLength": 100,
#                   "springConstant": 0.09,
#                   "avoidOverlap": 1
#                 },
#                 "minVelocity": 0.75,
#                 "solver": "forceAtlas2Based",
#                 "timestep": 0.22
#               }
#             }
#         '''
# 
#         net.set_options(opts)
#         # uncomment this to play with layout
#         #net.show_buttons(filter_=['physics'])
#         net.from_nx(g)
#         return net.show(name)
# =============================================================================

    
class Whole():
    def __init__(self, atomgroup, resolution=1, hyperres=False, 
                 neighbor_mask=np.ones((3,3,3))):
        self.voxels = Voxels(atomgroup, resolution, hyperres)
        self.voxels.label(neighbor_mask=neighbor_mask)
        self.bridges = Bridges(self.voxels.labels, self.voxels.nbox)
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
        box = dim2lattice(*shifted_atomgroup.dimensions)
        shifted_atomgroup.intersection(label_atomgroup).positions += shift @ box

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

    
class MDAWhole():
    """
    A MDAnalysis compatible Whole class which can be used
    as a just in time transformation or to make the complete
    trajectory whole and write the output.
    """
    def __init__(self, atomgroups, resolution=1):
        """
        Sets the atomgroup.
        """
        self.atomgroups = atomgroups
        self.resolution = resolution
    def __call__(self, ts):
        """
        Used in the MDA transformations logic for making a frame
        whole just in time.
        """
        for atomgroup in self.atomgroups:
            Whole(atomgroup, resolution=self.resolution)
        return ts
    @classmethod
    def whole_traj(cls, atomgroups, resolution=1, out='test.xtc', write_all=False, mol_whole=False):
        """
        Makes every frame whole, writes the xtc and returns the 
        whole atomgroup.
        """
        u = atomgroups[0].universe
        total_frames = len(u.trajectory)
        start_time = time()
        workflow = [MDAWhole(atomgroups, resolution)]
        if mol_whole:
            workflow.append(transformations.unwrap(u.atoms))
        u.trajectory.add_transformations(*workflow)
        combined_atomgroup = atomgroups[0]
        for atomgroup in atomgroups[1:]:
            combined_atomgroup += atomgroup
        with mda.Writer(out, combined_atomgroup.n_atoms) as W:
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
                if write_all:
                    W.write(u.atoms)
                else:
                    W.write(combined_atomgroup)
        print(f'\rDone, the whole thing took {(time()-start_time)/60:.2f} minutes.              ')
        return combined_atomgroup

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
    
    # OPIONAL
    required_grp.add_argument(
        '-x', '--trajectory', nargs='?', default=None, type=str,
        help='an MDAnalysis compatible coordinate or trajectory file (e.g. GRO, PDB, XTC)',
        )
    optional_grp.add_argument(
        '-sel', '--selections', nargs='?', default='not resname W WF ION', 
        type = str,
        help='the selection(s) to make whole, multiple selections can be separated with a semicolon (;)',
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
        '-wa', '--write_all', nargs='?', default=False, type=str,
        help='write all atoms from the original input (default=False)'
        )
    optional_grp.add_argument(
        '-mol', '--mol_whole', nargs='?', default=False, type=str,
        help='make molecules whole over pbc, this requires a tpr (default=False)'
        )
    optional_grp.add_argument(
    '-h', '--help', action="help",
    help='show this help message and exit',
    )
    # parse the arguments into the name space
    args = parser.parse_args()
    
    # Set the trajectory to the gro if it is not specified.
    if args.trajectory == None:
        args.out_file = 'whole.gro'
        args.trajectory = args.reference
    
    # Fixing the text input bools
    if args.write_all == False:
        pass
    elif args.write_all == "False":
        args.write_all = False
    elif args.write_all == "True":
        args.write_all = True
    else:
        raise ValueError('The -wa input should be either True or False.')
    if args.mol_whole == False:
        pass
    elif args.mol_whole == "False":
        args.mol_whole = False
    elif args.mol_whole == "True":
        args.mol_whole = True
    else:
        raise ValueError('The -mol input should be either True or False.')
    return args
    

def main():
    """
    Makes the trajectory whole and writes it.
    """
    args = read_arguments()
    selections = args.selections.split(';')
    print(selections)
    print(args)
    # Loading trajectory.
    if args.reference == args.trajectory:
        u = mda.Universe(args.reference, in_memory=False)    
    else:
        u = mda.Universe(args.reference, args.trajectory, in_memory=False)
    atomgroups = []
    for selection in selections:
        atomgroups.append(u.select_atoms(selection))
    # Make the complete trajectory whole and write it.
    MDAWhole.whole_traj(atomgroups, 
                        out=args.out_file, 
                        resolution=args.resolution,
                        write_all=args.write_all,
                        mol_whole=args.mol_whole)


if __name__ == "__main__":
    main()


# # Performance
# The main performance hit is the creation of the voxels2atoms dict. I could spent time on only checking half the neighbors due to symmetry, but it doesn't seem to be worth the effort. I think we could try no using a dict but the array itself, using np.unique and np.where(unique_value) instead. I am not sure if this will be faster in the end as you have to perform the np.where(unique_value) check for every unique value. This might turn out to be more expensive as the array gets bigger?
# 
# # To be added
# * A decent argparser.
# * Shifting non-selected particles with the shift of the label they are closest to (associative whole).
# * Submasking molecules e.g. segmenting only the tails, but projecting the shift on the whole lipid (sub whole). Sub whole could also work using the associate whole method, but I think it better is to make a separate function. Since it is cheaper if you do not have to do the distance check to know which label you belong to, as you are the same molecule!
