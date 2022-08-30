# mdvwhole
Density based object completion over all PBC. This repository will eventually be merged with [MDVoxelSegmentation](https://github.com/marrink-lab/MDVoxelSegmentation).

No copies of the box are made to complete the PBC, instead we use a graph based approach. Therefore object
completion is very light on memory and can cover an arbitrary amount of consecutive PBC crossings. The contact distance is based on
the voxel resolution (26 neighbors).

# Version
0.0.7.1 This version is safe to use for production. 

![alt text](https://user-images.githubusercontent.com/1488903/151573692-58d1eb6c-b6a2-444e-a7b8-937fa8ebc448.png)

## How to cite
For now it is best to cite the original MDVoxelSegmentation [article](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.1c00446) with the note that you used the new beta 'Whole' feature.

```
@article{Bruininks2021,
  doi = {10.1021/acs.jctc.1c00446},
  url = {https://doi.org/10.1021/acs.jctc.1c00446},
  year = {2021},
  month = oct,
  publisher = {American Chemical Society ({ACS})},
  volume = {17},
  number = {12},
  pages = {7873--7885},
  author = {Bart M. H. Bruininks and Albert S. Thie and Paulo C. T. Souza and Tsjerk A. Wassenaar and Shirin Faraji and Siewert J. Marrink},
  title = {Sequential Voxel-Based Leaflet Segmentation of Complex Lipid Morphologies},
  journal = {Journal of Chemical Theory and Computation}
}
```

## Install
Make sure you are installing for python 3.8 or newer.

`pip install mdvwhole`

Or clone this repository for the latest beta features and install using pip from inside the cloned folder.

`pip install -e .`

## Usage
Print help:

`mdvwhole -h`

Making a single gro whole:

`mdvwhole -f your_gro.gro`

Making a single gro whole and also complete single molecules:

`mdvwhole -f your_tpr.tpr -x your_gro.gro -mol True -o whole.gro`

Making a trajectory whole:

`mdvwhole -f your_gro.gro -x your_xtc.xtc -o whole.xtc`

Making a trajectory whole and also complete single molecules:

`mdvwhole -f your_tpr.tpr -x your_xtc.xtc -mol True -o whole.xtc`

Using a non-default selection:

`mdvwhole -f your_gro.gro -x your_xtc.xtc -o whole.xtc -sel 'not resname W WF ION and not name PO4'`

Using multiple selections to make whole:

`mdvwhole -f your_gro.gro -x your_xtc.xtc -o whole.xtc -sel 'resname POPC; not resname W WF ION POPC'`

Writing all atoms even if they were not included in a selection:

`mdvwhole -f your_gro.gro -x your_xtc.xtc -o whole.xtc -sel 'not resname W WF ION' -wa True`

Using an mdvoxelsegmentation 'clusters.npy' to make whole (can also be used without whole molecules and on a single gro).

`mdvwhole -f your_tpr.tpr -x your_xtc.xtc -o whole.xtc -wa True -mol True -clus your_clusters.npy`

https://user-images.githubusercontent.com/1488903/177659765-98287099-5619-4e45-b890-1de573437347.mp4

Using associative mdvwhole, a subselection of a molecules can be used for the selection. The displacement of the subselection is projected on the whole molecules:

`mdvwhole -f your_tpr.tpr -x your_xtc.xtc -o whole.xtc -sel 'name C3A C3B D3A D3B C4A C4B D4A D4B' -res 0.7 -wa True -mol True -asso True`

https://user-images.githubusercontent.com/1488903/177778628-ca61c694-fdd6-45f0-af78-644d63db9fe8.mp4

# Sneak peak at new feature, the aggregate shaped tile

https://user-images.githubusercontent.com/1488903/184648757-35735f05-bc7f-483e-b789-a74b834a1b6b.mp4

https://user-images.githubusercontent.com/1488903/185187491-0ba19d6e-cac3-460d-a97e-b1752060c050.mp4

https://user-images.githubusercontent.com/1488903/187438795-b5148b0d-f3de-4a4a-8f82-c45fbfb1b827.mp4



