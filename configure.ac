###########################################################################
# Prelude.
###########################################################################
AC_INIT([IBAMR],[080924],[griffith AT cims DOT nyu DOT edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config/IBAMR_config.h])
AM_INIT_AUTOMAKE([-Wall -Werror dist-bzip2 foreign])
AC_LANG(C++)

###########################################################################
# Checks for programs.
###########################################################################
AC_ARG_VAR(LIBS,[linker commands, e.g. -l<lib> if you require a nonstandard library <lib>])
AC_ARG_VAR(FCLIBS,[linker flags and commands, e.g. -L<lib dir> and -l<lib> for the "Fortran intrinsic and run-time libraries" that are required to successfully link a Fortran program or shared library])
AC_PROG_CC		# Check for a C compiler.
AC_PROG_CXX		# Check for a C++ compiler.
AC_PROG_FC		# Check for a Fortran compiler.
AC_FC_LIBRARY_LDFLAGS	# Determine the linker flags for the Fortran intrinsic and runtime libraries.
AC_FC_WRAPPERS
LIBS="$LIBS $FCLIBS"
AC_SUBST(F77, $FC)
AC_SUBST(FFLAGS, $FCFLAGS)
AC_SUBST(FLIBS, $FCLIBS)
AC_LIB_PREPARE_PREFIX
AC_LIB_RPATH

# Check for MPICC and MPICXX and configure MPI options.
AC_ARG_WITH([mpi],
  AS_HELP_STRING(--with-mpi=PATH,location of required MPI installation),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-mpi=PATH])
   fi
   if test -d "$withval" ; then
     PATH="${withval}:$PATH"
   fi
   if test -d "$withval/bin" ; then
     PATH="${withval}/bin:$PATH"
   fi],)
AC_ARG_WITH([mpicc],
  AS_HELP_STRING(--with-mpicc=ARG,MPI C compiler command),
  [MPICC=$withval],)
AC_ARG_WITH([mpicxx],
  AS_HELP_STRING(--with-mpicxx=ARG,MPI C++ compiler command),
  [MPICXX=$withval],)

AC_LANG_PUSH([C])
ACX_MPI(,AC_MSG_WARN([unable to find working MPI C complier wrapper]))
CC="$MPICC"
if test "$MPICC" != "" ; then
  CC="$MPICC"
fi
AC_LANG_POP([C])

ACX_MPI(,AC_MSG_WARN([unable to find working MPI C++ compiler wrapper]))
if test "$MPICXX" != "" ; then
  CXX="$MPICXX"
fi

if test "$MPILIBS" != "" ; then
  LIBS="$MPILIBS $LIBS"
fi

# Checks for C and C++ compiler features.
AC_HEADER_MAJOR
AC_HEADER_STDC

AC_CXX_NAMESPACES
if test "$ac_cv_cxx_namespaces" != "yes" ; then
  AC_MSG_ERROR([the specified C++ compiler does not appear to support namespaces])
fi

AC_CXX_HAVE_STL
if test "$ac_cv_cxx_have_stl" != "yes" ; then
  AC_MSG_ERROR([the specified C++ compiler does not appear support the Standard Template Library])
fi

# Checks for other programs.
AC_PROG_RANLIB

AC_ARG_WITH([M4],
  AS_HELP_STRING([--with-M4=ARG],[manually set M4 to ARG]))
AC_CHECK_PROGS(M4, $with_M4 gm4 m4,,)
if test -z "$M4" ; then
  AC_MSG_ERROR([m4 preprocessor not found, specify via --with-M4=ARG])
fi

AC_ARG_WITH([doxygen-dir],
  AS_HELP_STRING([--with-doxygen-dir=ARG],[manually set directory where the doxygen executable resides to ARG]),
  DOXYGEN_PATH=$withval
  AC_CHECK_PROGS(DOXYGEN, "$DOXYGEN_PATH/doxygen" doxygen,,),
  DOXYGEN_PATH=""
  AC_CHECK_PROGS(DOXYGEN, doxygen,,))
if test -z "$DOXYGEN" ; then
  AC_MSG_WARN([doxygen not found])
  echo "if doxygen is installed, specify its location via --with-doxygen-dir=ARG"
  HAVE_DOXYGEN="NO"
else
  HAVE_DOXYGEN="YES"
fi
AC_SUBST(HAVE_DOXYGEN,$HAVE_DOXYGEN)
AC_SUBST(DOXYGEN,$DOXYGEN)
AC_SUBST(DOXYGEN_PATH,$DOXYGEN_PATH)

AC_ARG_WITH([dot-dir],
  AS_HELP_STRING([--with-dot-dir=ARG],[manually set directory where the Graphviz dot executable resides to ARG]),
  DOT_PATH=$withval
  AC_CHECK_PROGS(DOT, "$DOT_PATH/dot" dot,,),
  DOT_PATH=""
  AC_CHECK_PROGS(DOT, dot,,))
if test -z "$DOT" ; then
  AC_MSG_WARN([dot not found])
  echo "if the Graphviz dot tool is installed, specify its location via --with-dot-dir=ARG"
  HAVE_DOT="NO"
else
  HAVE_DOT="YES"
fi
AC_SUBST(HAVE_DOT,$HAVE_DOT)
AC_SUBST(DOT,$DOT)
AC_SUBST(DOT_PATH,$DOT_PATH)

###########################################################################
# Checks for IBAMR configuration.
###########################################################################
echo
echo "========================="
echo "Configuring IBAMR options"
echo "========================="

# Checks to see whether large ghost cell widths should be enabled.
AC_ARG_ENABLE([large-ghost-cell-width],
  AS_HELP_STRING(--enable-large-ghost-cell-width,[enable the use of large ghost cell regions (this feature must be enabled for certain types of physical boundary condition routines to function properly) @<:@default=no@:>@]),
  [LARGE_GHOST_CELL_WIDTH=$enable_large_ghost_cell_width], [LARGE_GHOST_CELL_WIDTH=no])

if test "$LARGE_GHOST_CELL_WIDTH" == "yes" ; then
  echo "large ghost cell widths are ENABLED!"
  AC_DEFINE([USING_LARGE_GHOST_CELL_WIDTH],[true],[Boolean value indicating whether we are using large ghost cell widths])
else
  echo "large ghost cell widths are DISABLED!"
  echo "  certain physical boundary condition handling routines may not function properly"
  echo "  if you wish to enable support for large ghost cell widths, specify --enable-large-ghost-cell-width"
  AC_DEFINE([USING_LARGE_GHOST_CELL_WIDTH],[false],[Boolean value indicating whether we are using large ghost cell widths])
fi

# Checks to see whether expensive coarse-fine interface code should be
# enabled.
AC_ARG_ENABLE([expensive-cf-interpolation],
  AS_HELP_STRING(--enable-expensive-cf-interpolation,[enable the use of the expensive version of the quadratic coarse-fine interface interpolation code (this feature must be enabled for certain use cases, e.g. general anisotropic diffusion tensors) @<:@default=no@:>@]),
  [EXPENSIVE_CF_INTERPOLATION=$enable_expensive_cf_interpolation], [EXPENSIVE_CF_INTERPOLATION=no])

if test "$EXPENSIVE_CF_INTERPOLATION" == "yes" ; then
  echo "expensive version of the quadratic coarse-fine interface interpolation code is ENABLED!"
  AC_DEFINE([USING_EXPENSIVE_CF_INTERPOLATION],[true],[Boolean value indicating whether we are using the expensive version of the quadratic coarse-fine interface code])
else
  echo "expensive version of the quadratic coarse-fine interface interpolation code is DISABLED!"
  echo "  this feature must be enabled for certain use cases, e.g. general anisotropic diffusion tensors"
  echo "  if you wish to enable support for large ghost cell widths, specify --enable-expensive-cf-interpolation"
  AC_DEFINE([USING_EXPENSIVE_CF_INTERPOLATION],[false],[Boolean value indicating whether we are using the expensive version of the quadratic coarse-fine interface code])
fi

###########################################################################
# Checks for libraries.
###########################################################################

# Required external libraries.

# Check for BLAS and LAPACK libraries.
echo
echo "========================================"
echo "Configuring required package BLAS/LAPACK"
echo "========================================"
ACX_BLAS(,)
ACX_LAPACK(,)
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"

# Check for hdf5 library and configure hdf5 options.
echo
echo "================================="
echo "Configuring required package hdf5"
echo "================================="
AC_ARG_WITH([hdf5],
  AS_HELP_STRING(--with-hdf5=PATH,location of required hdf5 installation @<:@default=/usr/local/sfw/hdf5/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-hdf5=PATH])
   fi
   HDF5_DIR=$withval],
  [HDF5_DIR=/usr/local/sfw/hdf5/current])

if test -d "${HDF5_DIR}/lib" ; then
  LDFLAGS="-L${HDF5_DIR}/lib $LDFLAGS"
fi
if test -d "${HDF5_DIR}/include" ; then
  CPPFLAGS="-I${HDF5_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([hdf5.h],,AC_MSG_ERROR([could not find header file hdf5.h]))
AC_CHECK_HEADER([H5LT.h],,AC_MSG_ERROR([could not find header file H5LT.h]))

AC_LIB_HAVE_LINKFLAGS([hdf5],[z,sz])
LIBS="$LIBHDF5 $LIBS"

AC_LIB_HAVE_LINKFLAGS([hdf5_hl])
LIBS="$LIBHDF5_HL $LIBS"

# Check for hypre library and configure hypre options.
echo
echo "=================================="
echo "Configuring required package hypre"
echo "=================================="
AC_ARG_WITH([hypre],
  AS_HELP_STRING(--with-hypre=PATH,location of required hypre installation @<:@default=/usr/local/sfw/hypre/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-hypre=PATH])
   fi
   HYPRE_DIR=$withval],
  [HYPRE_DIR=/usr/local/sfw/hypre/current])

if test -d "${HYPRE_DIR}/lib" ; then
  LDFLAGS="-L${HYPRE_DIR}/lib $LDFLAGS"
fi
if test -d "${HYPRE_DIR}/include" ; then
  CPPFLAGS="-I${HYPRE_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([HYPRE.h],,AC_MSG_ERROR([could not find header file HYPRE.h]))

AC_LIB_HAVE_LINKFLAGS([HYPRE_utilities])
LIBS="$LIBHYPRE_UTILITIES $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_krylov])
LIBS="$LIBHYPRE_KRYLOV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_seq_mv])
LIBS="$LIBHYPRE_SEQ_MV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_parcsr_mv])
LIBS="$LIBHYPRE_PARCSR_MV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_parcsr_block_mv])
LIBS="$LIBHYPRE_PARCSR_BLOCK_MV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_IJ_mv])
LIBS="$LIBHYPRE_IJ_MV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_DistributedMatrix])
LIBS="$LIBHYPRE_DISTRIBUTEDMATRIX $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_MatrixMatrix])
LIBS="$LIBHYPRE_MATRIXMATRIX $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_Euclid])
LIBS="$LIBHYPRE_EUCLID $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_ParaSails])
LIBS="$LIBHYPRE_PARASAILS $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_DistributedMatrixPilutSolver])
LIBS="$LIBHYPRE_DISTRIBUTEDMATRIXPILUTSOLVER $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_parcsr_ls])
LIBS="$LIBHYPRE_PARCSR_LS $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_struct_mv])
LIBS="$LIBHYPRE_STRUCT_MV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_struct_ls])
LIBS="$LIBHYPRE_STRUCT_LS $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_sstruct_mv])
LIBS="$LIBHYPRE_SSTRUCT_MV $LIBS"

AC_LIB_HAVE_LINKFLAGS([HYPRE_sstruct_ls])
LIBS="$LIBHYPRE_SSTRUCT_LS $LIBS"

# Check for PETSc library and configure PETSc options.
echo
echo "=================================="
echo "Configuring required package PETSc"
echo "=================================="
AC_ARG_VAR(PETSC_DIR,[the location of the PETSc installation that is to be used.  Note that multiple versions of PETSc may be installed within the same filesystem, with each version corresponding to a different value of PETSC_DIR.])
AC_ARG_VAR(PETSC_ARCH,[the PETSc configuration that is to be used, corresponding to the configuration located in ${PETSC_DIR}/${PETSC_ARCH}.  Note that multiple configurations of PETSc can be installed for a particular version of the PETSc library.  Each configuration will correspond to a different value of PETSC_ARCH.])
AC_ARG_WITH([petsc],
  AS_HELP_STRING(--with-petsc=PATH,location of required PETSc installation @<:@default=PETSC_DIR@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-petsc=PATH])
   fi
   PETSC_DIR=$withval],)
AC_ARG_WITH([petsc-arch],
  AS_HELP_STRING(--with-petsc-arch=ARCH,PETSc architecture @<:@default=PETSC_ARCH@:>@),
  [PETSC_ARCH=$withval],)
AC_SUBST(PETSC_DIR,$PETSC_DIR)
AC_SUBST(PETSC_ARCH,$PETSC_ARCH)

echo "using PETSc architecture = ${PETSC_ARCH}"
echo "using PETSc root    directory = ${PETSC_DIR}"
echo "using PETSc lib     directory = ${PETSC_DIR}/${PETSC_ARCH}/lib"
echo "using PETSc include directory = ${PETSC_DIR}/include, ${PETSC_DIR}/${PETSC_ARCH}/include"

if test -d "${PETSC_DIR}/${PETSC_ARCH}/lib" ; then
  LDFLAGS="-L${PETSC_DIR}/${PETSC_ARCH}/lib $LDFLAGS"
fi
if test -d "${PETSC_DIR}/${PETSC_ARCH}/include" ; then
  CPPFLAGS="-I${PETSC_DIR}/${PETSC_ARCH}/include $CPPFLAGS"
fi
if test -d "${PETSC_DIR}/include" ; then
  CPPFLAGS="-I${PETSC_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([petsc.h],,AC_MSG_ERROR([could not find header file petsc.h]))

AC_LIB_HAVE_LINKFLAGS([petsc])
LIBS="$LIBPETSC $LIBS"

AC_LIB_HAVE_LINKFLAGS([petscvec])
LIBS="$LIBPETSCVEC $LIBS"

AC_LIB_HAVE_LINKFLAGS([petscmat])
LIBS="$LIBPETSCMAT $LIBS"

AC_LIB_HAVE_LINKFLAGS([petscdm])
LIBS="$LIBPETSCDM $LIBS"

AC_LIB_HAVE_LINKFLAGS([petscksp])
LIBS="$LIBPETSCKSP $LIBS"

AC_LIB_HAVE_LINKFLAGS([petscsnes])
LIBS="$LIBPETSCSNES $LIBS"

AC_LIB_HAVE_LINKFLAGS([petscts])
LIBS="$LIBPETSCTS $LIBS"

# Check for SAMRAI library and configure SAMRAI options.
echo
echo "==================================="
echo "Configuring required package SAMRAI"
echo "==================================="
AC_ARG_WITH([samrai],
  AS_HELP_STRING(--with-samrai=PATH,location of required SAMRAI installation @<:@default=/usr/local/sfw/samrai/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-samrai=PATH])
   fi
   SAMRAI_DIR=$withval],
  [SAMRAI_DIR=/usr/local/sfw/samrai/current])

if test -d "${SAMRAI_DIR}/lib" ; then
  LDFLAGS="-L${SAMRAI_DIR}/lib $LDFLAGS"
fi
if test -d "${SAMRAI_DIR}/include" ; then
  CPPFLAGS="-I${SAMRAI_DIR}/include $CPPFLAGS"
fi

AC_SUBST(SAMRAI_DIR,[${SAMRAI_DIR}])
AC_SUBST(SAMRAI_FORTDIR,[${SAMRAI_DIR}/include])

AC_CHECK_HEADER([SAMRAI_config.h],,AC_MSG_ERROR([could not find header file SAMRAI_config.h]))

AC_LIB_HAVE_LINKFLAGS([SAMRAI])
LIBS="$LIBSAMRAI $LIBS"

AC_ARG_ENABLE([samrai-2d],
  AS_HELP_STRING(--enable-samrai-2d,enable optional support for two-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI2D=$enablevar], [SAMRAI2D=yes])
AM_CONDITIONAL([SAMRAI2D_ENABLED],[test "$SAMRAI2D" == "yes"])

if test "$SAMRAI2D" == "yes"; then
  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_hier])
  LIBS="$LIBSAMRAI2D_HIER $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_xfer])
  LIBS="$LIBSAMRAI2D_XFER $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_pdat_std])
  LIBS="$LIBSAMRAI2D_PDAT_STD $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_math_std])
  LIBS="$LIBSAMRAI2D_MATH_STD $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_mesh])
  LIBS="$LIBSAMRAI2D_MESH $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_geom])
  LIBS="$LIBSAMRAI2D_GEOM $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_solv])
  LIBS="$LIBSAMRAI2D_SOLV $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_algs])
  LIBS="$LIBSAMRAI2D_ALGS $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI2d_appu])
  LIBS="$LIBSAMRAI2D_APPU $LIBS"
else
  AC_MSG_WARN([configuring without the two-dimensional SAMRAI library])
fi

AC_ARG_ENABLE([samrai-3d],
  AS_HELP_STRING(--enable-samrai-3d,enable optional support for three-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI3D=$enablevar], [SAMRAI3D=yes])
AM_CONDITIONAL([SAMRAI3D_ENABLED],[test "$SAMRAI3D" == "yes"])

if test "$SAMRAI3D" == "yes"; then
  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_hier])
  LIBS="$LIBSAMRAI3D_HIER $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_xfer])
  LIBS="$LIBSAMRAI3D_XFER $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_pdat_std])
  LIBS="$LIBSAMRAI3D_PDAT_STD $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_math_std])
  LIBS="$LIBSAMRAI3D_MATH_STD $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_mesh])
  LIBS="$LIBSAMRAI3D_MESH $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_geom])
  LIBS="$LIBSAMRAI3D_GEOM $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_solv])
  LIBS="$LIBSAMRAI3D_SOLV $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_algs])
  LIBS="$LIBSAMRAI3D_ALGS $LIBS"

  AC_LIB_HAVE_LINKFLAGS([SAMRAI3d_appu])
  LIBS="$LIBSAMRAI3D_APPU $LIBS"
else
  AC_MSG_WARN([configuring without the three-dimensional SAMRAI library])
fi

# Optional external libraries.

# Check for Blitz++ library and configure Blitz++ options.
AC_ARG_ENABLE([blitz],
  AS_HELP_STRING(--enable-blitz,enable support for the optional Blitz++ library @<:@default=yes@:>@),
  [BLITZ_ENABLED=$enableval], [BLITZ_ENABLED=yes])
if test "$BLITZ_ENABLED" == "" ; then
  BLITZ_ENABLED=yes
fi

AC_ARG_WITH([blitz],
  AS_HELP_STRING(--with-blitz=PATH,location of optional Blitz++ installation @<:@default=/usr/local/sfw/blitz/current@:>@),
  [if test "$BLITZ_ENABLED" == "no" ; then
     AC_MSG_WARN([--with-blitz is specified, but support for Blitz++ is DISABLED])
     BLITZ_DIR=NONE
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([you must specify an existing directory when using --with-blitz=PATH])
     fi
     BLITZ_DIR=$withval
   fi],
  [BLITZ_DIR=/usr/local/sfw/blitz/current])

if test "$BLITZ_ENABLED" == "yes"; then
  echo
  echo "===================================="
  echo "Configuring optional package Blitz++"
  echo "===================================="

  if test -d "${BLITZ_DIR}/lib" ; then
    LDFLAGS="-L${BLITZ_DIR}/lib $LDFLAGS"
  fi
  if test -d "${BLITZ_DIR}/include" ; then
   CPPFLAGS="-I${BLITZ_DIR}/include $CPPFLAGS"
  fi

  AC_CHECK_HEADER([blitz/blitz.h],,AC_MSG_WARN([could not find header file blitz/blitz.h]))
  AC_LIB_HAVE_LINKFLAGS([blitz])
  LIBS="$LIBBLITZ $LIBS"
else
  echo
  echo "===================================="
  echo "Optional package Blitz++ is DISABLED"
  echo "===================================="
fi

# Check for Silo library and configure Silo options.
AC_ARG_ENABLE([silo],
  AS_HELP_STRING(--enable-silo,enable support for the optional Silo library @<:@default=yes@:>@),
  [SILO_ENABLED=$enableval], [SILO_ENABLED=yes])
if test "$SILO_ENABLED" == "" ; then
  SILO_ENABLED=yes
fi

AC_ARG_WITH([silo],
  AS_HELP_STRING(--with-silo=PATH,location of optional Silo installation @<:@default=/usr/local/sfw/silo/current@:>@),
  [if test "$SILO_ENABLED" == "no" ; then
     AC_MSG_WARN([--with-silo is specified, but support for Silo is DISABLED])
     SILO_DIR=NONE
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([you must specify an existing directory when using --with-silo=PATH])
     fi
     SILO_DIR=$withval
   fi],
  [SILO_DIR=/usr/local/sfw/silo/current])

if test "$SILO_ENABLED" == "yes"; then
  echo
  echo "================================="
  echo "Configuring optional package Silo"
  echo "================================="

  if test -d "${SILO_DIR}/lib" ; then
    LDFLAGS="-L${SILO_DIR}/lib $LDFLAGS"
  fi
  if test -d "${SILO_DIR}/include" ; then
    CPPFLAGS="-I${SILO_DIR}/include $CPPFLAGS"
  fi

  AC_CHECK_HEADER([silo.h],,AC_MSG_WARN([could not find header file silo.h]))
  AC_LIB_HAVE_LINKFLAGS([silo])
  if test "$HAVE_LIBSILO" == "no" ; then
    AC_MSG_WARN([could not find libsilo; checking for libsiloh5])
    AC_LIB_HAVE_LINKFLAGS([siloh5])
    LIBS="$LIBSILOH5 $LIBS"
  fi
else
  echo
  echo "================================="
  echo "Optional package Silo is DISABLED"
  echo "================================="
fi

###########################################################################
# Output files.
###########################################################################
echo
echo "================"
echo "Outputting files"
echo "================"
AC_CONFIG_SUBDIRS([ibtk])
AC_CONFIG_FILES([
  Makefile
  config/make.inc
  doc/Makefile
  doc/Doxyfile
  sample/IB/Makefile
  sample/empty_project/Makefile
  src/Makefile
  src/fortran/Makefile
  src/IB/Makefile
  src/adv_diff/Makefile
  src/adv_diff/fortran/Makefile
  src/advect/Makefile
  src/advect/fortran/Makefile
  src/navier_stokes/Makefile
  src/navier_stokes/fortran/Makefile
  src/utilities/Makefile
  src/utilities/get_version.C
  src/tests/Makefile
  src/tests/IB/Makefile
  src/tests/IB_multilevel/Makefile
  src/tests/IB_staggered/Makefile
  src/tests/IB_target/Makefile
  src/tests/IB_target_staggered/Makefile
  src/tests/adv_diff/Makefile
  src/tests/advect/Makefile
  src/tests/navier_stokes/Makefile
  src/tests/navier_stokes/analytic_soln/Makefile
  src/tests/navier_stokes/analytic_soln_staggered/Makefile
  src/tests/navier_stokes/lid_driven_cavity/Makefile
  src/tests/navier_stokes/lid_driven_cavity_staggered/Makefile
  src/tests/navier_stokes/periodic_shear_layer/Makefile
  src/tests/navier_stokes/periodic_shear_layer_staggered/Makefile
  src/tests/navier_stokes/regularized_lid_driven_cavity_staggered/Makefile
  src/tests/pIB/Makefile
  src/tools/Makefile
])
AC_OUTPUT
