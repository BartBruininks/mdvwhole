# mdvwhole
Density based object completion over PBC. This repository eventually will be merged with MDVoxelSegmentation.


## Usage
Install the dependencies at the top of whole.py using pip.

Making a single gro whole:

`python whole.py -f your_gro.gro -x your_gro -o whole.gro`

Making a trajectory whole:

`python whole.py -f your_gro.gro -x your_xtc -o whole.xtc`

Using a non-default selection:

`python whole.py -f your_gro.gro -x your_xtc -o whole.xtc -sel 'not resname W WF ION and not name PO4'`
