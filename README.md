# Poiseuille

The 2D adaptation is done with hanging nodes, till now, it is just possible to do 1  adapt step. The reason is only adapt if triangles is done, so for the new polygons,  the algorithm has not clear how to adapt them.

On the other hand,  first adapt in 2D seems to have problems, since values of solution are inf, and ctk files re not read by paraview.
