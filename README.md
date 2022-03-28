# mdvwhole
Density based object completion over orthorhombic PBC. All PBC will be supported in the future, but for now it is not! This repository will eventually be merged with [MDVoxelSegmentation](https://github.com/marrink-lab/MDVoxelSegmentation).

No copies of the box are made to complete the PBC, instead we use a graph based approach. Therefore object
completion is very light on memory and can cover an arbitrary amount of consecutive PBC crossings. The contact distance is based on
the voxel resolution (26 neighbors).

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

If you are running into issues with the astype('int32') statement, it means you are probably running in a lower version of python, or you do not have te latest version of numba. By commenting out the @jit lines one can disable the numba acceleration. This will make it roughly 1.5x slower, but now it will run on any version of python3. It would be nice if this was a flag which could be set.

## Usage
Making a single gro whole:

`mdvwhole -f your_gro.gro -x your_gro -o whole.gro`

Making a trajectory whole:

`mdvwhole -f your_gro.gro -x your_xtc -o whole.xtc`

Using a non-default selection:

`mdvwhole -f your_gro.gro -x your_xtc -o whole.xtc -sel 'not resname W WF ION and not name PO4'`
