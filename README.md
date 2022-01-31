# mdvwhole
Density based object completion over PBC. This repository will eventually be merged with [MDVoxelSegmentation](https://github.com/marrink-lab/MDVoxelSegmentation).

No copies of the box are made to complete the PBC, instead we use a graph based approach. Therefore object
completion is very light on memory and can cover an arbitrary amount of consecutive PBC crossings.

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
`pip install mdvwhole`

## Usage
Making a single gro whole:

`mdvwhole -f your_gro.gro -x your_gro -o whole.gro`

Making a trajectory whole:

`mdvwhole -f your_gro.gro -x your_xtc -o whole.xtc`

Using a non-default selection:

`mdvwhole -f your_gro.gro -x your_xtc -o whole.xtc -sel 'not resname W WF ION and not name PO4'`

![alt text](https://user-images.githubusercontent.com/1488903/151573692-58d1eb6c-b6a2-444e-a7b8-937fa8ebc448.png)
