This example demonstrates flow past a cylinder using the IBFE method. The cylinder can be modeled in two ways:

1. Using only a surface mesh (set `USE_BOUNDARY_MESH=TRUE` in the input file)
2. Using a mesh that approximates both the boundary and interior of the cylinder (set `USE_BOUNDARY_MESH=FALSE`)

No-slip boundary conditions on the cylinder are enforced using penalty springs and/or a frictional penalty. The drag and lift coefficients of the cylinder are computed and output at intervals specified by `data_dump_interval`.

For full details, see: Hybrid finite difference/finite element immersed boundary method. *International Journal for Numerical Methods in Biomedical Engineering*, 33(11):e2888, 2017.
