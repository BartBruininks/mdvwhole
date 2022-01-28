# mdvwhole
Density based object completion over PBC. This repository eventually will be merged with MDVoxelSegmentation.

## Install
pip install mdvwhole

## Usage
Install the dependencies at the top of whole.py using pip.

Making a single gro whole:

`mdvwhole -f your_gro.gro -x your_gro -o whole.gro`

Making a trajectory whole:

`mdvwhole -f your_gro.gro -x your_xtc -o whole.xtc`

Using a non-default selection:

`mdvwhole -f your_gro.gro -x your_xtc -o whole.xtc -sel 'not resname W WF ION and not name PO4'`
