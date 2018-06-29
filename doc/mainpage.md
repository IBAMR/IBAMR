IBAMR
=====

An adaptive and distributed-memory parallel implementation of the immersed boundary (IB) method

What Is %IBAMR?
--------------

%IBAMR is a distributed-memory parallel implementation of the immersed boundary (IB) method with support for Cartesian grid adaptive mesh refinement (AMR).  Support for distributed-memory parallelism is via <a href="http://www.mcs.anl.gov/research/projects/mpi" target="_blank">MPI</a>, the Message Passing Interface.

Core %IBAMR functionality relies upon several high-quality open-source libraries, including:
 * <a href="https://computation-rnd.llnl.gov/SAMRAI" target="_blank">SAMRAI</a>, the Structured Adaptive Mesh Refinement Application Infrastructure
 * <a href="http://www.mcs.anl.gov/petsc" target="_blank">PETSC</a>, the Portable, Extensible Toolkit for Scientific Computation
 * <a href="http://libmesh.sourceforge.net" target="_blank">libMesh</a>, a C++ finite element library
 * <a href="http://computation.llnl.gov/casc/linear_solvers/sls_hypre.html" target="_blank">hypre</a>, a library of high performance preconditioners that features parallel multigrid methods for both structured and unstructured grid problems.

%IBAMR also uses functionality provided by a number of additional third-party libraries, including: <a href="http://www.boost.org" target="_blank">Boost</a>; <a href="http://eigen.tuxfamily.org/index.php" target="_blank">Eigen</a>; <a href="http://www.hdfgroup.org/HDF5" target="_blank">HDF5</a>; <a href="http://muparser.beltoforion.de" target="_blank">muParser</a>; and <a href="https://wci.llnl.gov/codes/silo" target="_blank">Silo</a>.

%IBAMR outputs visualization files that can be read by the <a href="https://wci.llnl.gov/codes/visit" target="_blank">VisIt Visualization Tool</a>.

What Is the IB Method?
----------------------

The immersed boundary (IB) method is a general-purpose numerical method for simulating fluid-structure interaction.  The IB formulation of such problems uses an Eulerian description of the fluid and a Lagrangian description of the structure.  Interaction equations that couple the Eulerian and Lagrangian variables take the form of integral equations with delta function kernels.

For general information about the IB method, see <a href="http://math.nyu.edu/faculty/peskin" target="_blank"></a>.  For visualizations of simulations that use %IBAMR, see <a href="http://cims.nyu.edu/~griffith" target="_blank"></a>.


Support
-------

Support for %IBAMR is available via the <a href="https://groups.google.com/forum/#!forum/ibamr-users" target="_blank">IBAMR Users Google Group</a> (ibamr-users@googlegroups.com).  Discussion related to the continued development of %IBAMR is via the <a href="http://groups.google.com/group/ibamr-dev" target="_blank">IBAMR Developers Google Group</a> (ibamr-dev@googlegroups.com). 

Bugs and Other Issues
---------------------

Please use the GitHub issue tracking system to report bugs, feature requests, or
 other issues with %IBAMR.

Acknowledgments
---------------

%IBAMR development is supported in part by an NSF <i>Software Infrastructure for Sustained Innovation</i> award (NSF OCI 1047734).  Work to extend %IBAMR to support finite element mechanics models is also supported in part by the NSF (NSF DMS 1016554).  We gratefully acknowledge this support.

