# mdvwhole
Density based object completion over PBC. This repository will eventually be merged with MDVoxelSegmentation.

No copies of the box are made to complete the PBC, instead we use a graph base approach. Therefore object
completion is very light on memory and can cover an arbitrary amount of consecutive PBC crossings.

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
