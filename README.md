# WARNING
We found a functionality breaking bug for small molecules crossing boundaries at diagonals. We are looking into this. For now, use MDVWhole with care, if it works for you, good. But we know the output is wrong in some cases. We will fix this ASAP and apologize for any inconvenience this might have brought on our users.

# mdvwhole
Density based object completion over all PBC, the [manuscript](https://doi.org/10.1021/acs.jcim.2c01574) is published in JCIM.

No copies of the box are made to complete the PBC, instead we use a graph based approach. Therefore object
completion is very light on memory and can cover an arbitrary amount of consecutive PBC crossings. The contact distance is based on
the voxel resolution (26 neighbors). This repository is closely related to [MDVoxelSegmentation](https://github.com/marrink-lab/MDVoxelSegmentation) and in the future I hope to merge them.

# Version
0.0.7.1 This version is safe to use for production. 

![alt text](https://user-images.githubusercontent.com/1488903/151573692-58d1eb6c-b6a2-444e-a7b8-937fa8ebc448.png)

## How to cite
```
@article{Bruininks2023,
  doi = {10.1021/acs.jcim.2c01574},
  url = {https://doi.org/10.1021/acs.jcim.2c01574},
  year = {2023},
  month = may,
  publisher = {American Chemical Society ({ACS})},
  author = {Bart M. H. Bruininks and Tsjerk A. Wassenaar and Ilpo Vattulainen},
  title = {Unbreaking Assemblies in Molecular Simulations with Periodic Boundaries},
  journal = {Journal of Chemical Information and Modeling}
}
```

## Install
Make sure you are installing for python 3.8 or newer.

`pip install git+https://github.com/BartBruininks/mdvwhole`

Or clone this repository for the latest beta features and install using pip from inside the cloned folder.

`pip install -e .`

## Usage
Print help:

`mdvwhole -h`

Making a single gro whole:

`mdvwhole -f your_gro.gro`

Using a custom resolution (setting the resolution to 0.8 nm and ignore the voxel scaling artifact warning by using a minus sign in fron the the resolution):

`mdvwhole -f your_gro.gro -x your_xtc.xtc -res -0.4` # 0.4 a lower value than this leads to gaps for CG Martini models

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

Using associative mdvwhole, a subselection of a molecules can be used for the selection. The displacement of the subselection is projected on the whole molecules. A negative sign in front of the resolution allows for large deviations from the set resolution if required. Otherwise a maximum deviation of 5% is tolerated:

`mdvwhole -f your_tpr.tpr -x your_xtc.xtc -o whole.xtc -sel 'name C3A C3B D3A D3B C4A C4B D4A D4B' -res -0.7 -wa True -mol True -asso True`

https://user-images.githubusercontent.com/1488903/177778628-ca61c694-fdd6-45f0-af78-644d63db9fe8.mp4

If an object starts to get self contact over the PBC (end of the video) it can no longer be completed with MDVWhole. We see this as a feature, for such unwanted contacts are most often a strong indication that the box used should be bigger.

https://user-images.githubusercontent.com/1488903/203772841-8f70100a-51f0-418a-8e2c-c4aa86ed1284.mp4

# Sneak peak at new feature, the aggregate shaped tile

https://user-images.githubusercontent.com/1488903/184648757-35735f05-bc7f-483e-b789-a74b834a1b6b.mp4

https://user-images.githubusercontent.com/1488903/185187491-0ba19d6e-cac3-460d-a97e-b1752060c050.mp4

https://user-images.githubusercontent.com/1488903/187438795-b5148b0d-f3de-4a4a-8f82-c45fbfb1b827.mp4



