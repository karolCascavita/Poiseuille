# Poiseuille

The 2D adaptation is done with hanging nodes, till now, it is just possible to do 1  adaptative step. The reason is only adaptation is performe for triangles, so for the new polygons the algorithm has not clear how to adapt them.

On the other hand,  first adaptation in 2D seems to have problems, since values of solution are inf, and vtk files re not read by paraview.
