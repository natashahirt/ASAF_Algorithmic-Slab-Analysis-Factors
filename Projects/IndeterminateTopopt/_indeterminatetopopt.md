# Indeterminate Analysis and Topology Optimization

This project was performed in April-May 2025, and is featured in the final chapter of my SMBT thesis. It consists of two parts:

1. an alternative analysis method to discrete tributary area + Hillerborg strip analysis that takes into consideration the indeterminate nature of concrete floor systems
2. a "topology optimization" loop that outputs mediocre results given a starting ground structure

One of the advantages of the indeterminate analysis, which distributes point loads (i.e. from mesh divisions) to each element in the layout that has a valid load path to the point load (e.g. in the case of isotropic slabs, an element has a valid load path if a line drawn normal to the beam intersects with the point load at some point). The percentage of the total point load applied to each beam is proportional to the stiffness of the beam and the distance from the beam to the point load relative to other load paths available to that point load.

The main to-do item in this project is to implement directionality in the indeterminate load assignment (i.e. unidirectional and bidirectional.)
