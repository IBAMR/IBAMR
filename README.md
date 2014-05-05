IBAMR
=====

An adaptive and distributed-memory parallel implementation of the immersed boundary (IB) method

What Is IBAMR?
--------------

IBAMR is a distributed-memory parallel implementation of the immersed boundary (IB) method with support for Cartesian grid adaptive mesh refinement (AMR).  Support for distributed-memory parallelism is via [MPI](http://www.mcs.anl.gov/research/projects/mpi), the Message Passing Interface.

Core IBAMR functionality relies upon several high-quality open-source libraries, including:
 * [SAMRAI](https://computation-rnd.llnl.gov/SAMRAI), the Structured Adaptive Mesh Refinement Application Infrastructure
 * [PETSc](http://www.mcs.anl.gov/petsc), the Portable, Extensible Toolkit for Scientific Computation
 * [libMesh](http://libmesh.sourceforge.net), a C++ finite element library
 * [*hypre*](http://computation.llnl.gov/casc/linear_solvers/sls_hypre.html), a library of high performance preconditioners that features parallel multigrid methods for both structured and unstructured grid problems.
 
IBAMR also uses functionality provided by a number of additional third-party libraries, including: [Boost](http://www.boost.org); [Eigen](http://eigen.tuxfamily.org/index.php); [HDF5](http://www.hdfgroup.org/HDF5); [muParser](http://muparser.beltoforion.de); and [Silo](https://wci.llnl.gov/codes/silo).

IBAMR outputs visualization files that can be read by the [VisIt Visualization Tool](https://wci.llnl.gov/codes/visit).

What Is the IB Method?
----------------------

The immersed boundary (IB) method is a general-purpose numerical method for simulating fluid-structure interaction.  The IB formulation of such problems uses an Eulerian description of the fluid and a Lagrangian description of the structure.  Interaction equations that couple the Eulerian and Lagrangian variables take the form of integral equations with delta function kernels. </p><p>For general information about the IB method, see [http://math.nyu.edu/faculty/peskin].  For visualizations of simulations that use IBAMR, see [http://cims.nyu.edu/~griffith].

Getting Started
---------------

1. [Download and build the required third-party libraries](/p/ibamr/wiki/Building_Third_Party_Libraries)
2. [Download and build IBAMR](/p/ibamr/wiki/Building_IBAMR)

Source code documentation [for IBAMR is available on-line](http://ibamr.googlecode.com/svn/doc/ibamr/HEAD/html/annotated.html).  Source code documentation is also available [for the IBTK support library](http://ibamr.googlecode.com/svn/doc/ibtk/HEAD/html/annotated.html).  File format documentation [is also available on-line](http://ibamr.googlecode.com/svn/doc/ibamr/HEAD/html/classIBAMR_1_1IBStandardInitializer.html#_details).

There are some basic instructions on using or building IBAMR at NYU:
  * [Using IBAMR on the CIMS Linux network](/p/ibamr/wiki/IBAMR_CIMS)
  * [Building IBAMR on the CIMS Linux network](/p/ibamr/wiki/IBAMR_CIMS_Build)
  * [Building IBAMR on the NYU cardiac HPC cluster](/p/ibamr/wiki/IBAMR_cardiac_Build)
  * [Building IBAMR on the NYU bowery HPC cluster](/p/ibamr/wiki/IBAMR_bowery_Build)

There are also some basic instructions on building IBAMR on various platforms:
  * [Building IBAMR on a system running Linux](/p/ibamr/wiki/IBAMR_Linux_Build)
  * [Building IBAMR on a system running Mac OS X](/p/ibamr/wiki/IBAMR_OS_X_Build)

There are also [guidelines on how to contribute to the IBAMR project](/p/ibamr/wiki/IBAMR_Development).

Support
-------

Support for IBAMR is available via the [IBAMR Users](http://groups.google.com/group/ibamr-users) (ibamr-users) Google Group.  Discussion related to the continued development of IBAMR is via the [IBAMR Developers](http://groups.google.com/group/ibamr-dev) (ibamr-dev) Google Group. There is also list of [frequently asked questions](/p/ibamr/wiki/FAQ).

Bugs and Other Issues
---------------------

Please use the GitHub issue tracking system to report bugs, feature requests, or other issues with IBAMR.

Acknowledgments
---------------

IBAMR development is supported in part by an NSF <i>Software Infrastructure for Sustained Innovation</i> award (NSF OCI 1047734).  Work to extend IBAMR to support finite element mechanics models is also supported in part by the NSF (NSF DMS 1016554).  We gratefully acknowledge this support.
