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

The immersed boundary (IB) method is a general-purpose numerical method for simulating fluid-structure interaction.  The IB formulation of such problems uses an Eulerian description of the fluid and a Lagrangian description of the structure.  Interaction equations that couple the Eulerian and Lagrangian variables take the form of integral equations with delta function kernels.

For general information about the IB method, see [math.nyu.edu/faculty/peskin](http://math.nyu.edu/faculty/peskin).  For additional information about the IBAMR software, see [ibamr.github.io](http://ibamr.github.io).  We are happy to host visualizations of simulations that use IBAMR.

Getting Started
---------------

IBAMR requires a number of [third-party libraries](../../wiki/ThirdPartyLibraries).  [Sample build instructions are provided](../../wiki/Building) for a typical Linux installation.

Documentation
-------------

Source code documentation [for the IBAMR library and supporting IBTK library is available on-line](http://ibamr.github.io/docs). File format documentation [is also available on-line](https://ibamr.github.io/IBAMR-docs/ibamr/html/class_i_b_a_m_r_1_1_i_b_standard_initializer.html#details).  There is also list of [frequently asked questions](https://ibamr.github.io/faq).

Support
-------

Support for IBAMR is available via the [IBAMR Users Google Group (ibamr-users@googlegroups.com)](http://groups.google.com/group/ibamr-users).  Discussion related to the continued development of IBAMR is via the [IBAMR Developers Google Group (ibamr-dev@googlegroups.com)](http://groups.google.com/group/ibamr-dev).

Bugs and Other Issues
---------------------

Please use the GitHub issue tracking system to report bugs, feature requests, or other issues with IBAMR.

Acknowledgments
---------------

IBAMR development is supported in part by NSF <i>Software Infrastructure for Sustained Innovation</i> awards OAC 1450327 (to UNC-Chapel Hill), OAC 1450374 (to Northwestern University), and OAC 1607042 (to Rice University).  Additional support is provided by NSF CAREER award OAC 1652541 (to UNC-Chapel Hill).  Prior support was provided by NSF awards DMS 1016554 (to New York University), DMS 1460368 (to UNC-Chapel Hill), OAC 1047734 (to New York University), and OAC 1460334 (to UNC-Chapel Hill).  We gratefully acknowledge this support.
