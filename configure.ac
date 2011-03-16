###########################################################################
# Prelude.
###########################################################################
AC_PREREQ([2.68])
AC_INIT([IBAMR],[101029],[griffith AT cims DOT nyu DOT edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config/IBAMR_config.h])
AM_INIT_AUTOMAKE([1.11.1 -Wall -Werror dist-bzip2 foreign])
AM_MAINTAINER_MODE([disable])
AC_LANG(C++)

###########################################################################
# Checks for programs.
###########################################################################
AC_ARG_VAR(LDFLAGS,[linker flags, e.g. -L<lib dir> if you have libraries in a nonstandard directory <lib dir>])
AC_ARG_VAR(LIBS,[linker commands, e.g. -l<lib> if you require a nonstandard library <lib>])
AC_ARG_VAR(FCLIBS,[linker flags and commands, e.g. -L<lib dir> and -l<lib> for the "Fortran intrinsic and run-time libraries" that are required to successfully link a Fortran program or shared library])
AC_PROG_CC		# Check for a C compiler.
AC_PROG_CXX		# Check for a C++ compiler.
AC_PROG_CXXCPP          # Check for a C++ preprocessor.
AC_PROG_FC		# Check for a Fortran compiler.
AC_PROG_F77		# Check for a Fortran 77 compiler.
AC_FC_LIBRARY_LDFLAGS	# Determine the linker flags for the Fortran intrinsic and runtime libraries.
AC_FC_WRAPPERS
LIBS="$LIBS $FCLIBS"

if test "x$F77" != "x" -a "$F77" != "$FC" ; then
  AC_MSG_WARN([environment variables F77 and FC differ; using F77=$FC])
fi
F77=$FC

if test "x$FFLAGS" != "x" -a "$FFLAGS" != "$FCFLAGS" ; then
  AC_MSG_WARN([environment variables FFLAGS and FCFLAGS differ; using FFLAGS=$FCFLAGS])
fi
FFLAGS=$FCFLAGS

if test "x$FLIBS" != "x" -a "$FLIBS" != "$FCLIBS" ; then
  AC_MSG_WARN([environment variables FCLIBS and FLIBS differ; using FLIBS=$FCLIBS])
fi
FLIBS=$FCLIBS

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
   fi
   if test -d "$withval/include" ; then
     CPPFLAGS="-I${withval}/include $CPPFLAGS"
   fi],)
AC_ARG_WITH([mpicc],
  AS_HELP_STRING(--with-mpicc=ARG,MPI C compiler command),
  [MPICC=$withval],)
AC_ARG_WITH([mpicxx],
  AS_HELP_STRING(--with-mpicxx=ARG,MPI C++ compiler command),
  [MPICXX=$withval],)

AC_LANG_PUSH([C])
AX_MPI([],[AC_MSG_WARN([unable to find working MPI C complier wrapper])])
AC_LANG_POP([C])

AX_MPI([],[AC_MSG_WARN([unable to find working MPI C++ compiler wrapper])])

if test "$MPILIBS" != "" ; then
  LIBS="$MPILIBS $LIBS"
fi

# Checks for C and C++ compiler features.
AC_HEADER_MAJOR
AC_HEADER_STDC

AX_CXX_NAMESPACES
if test "$ax_cv_cxx_namespaces" != "yes" ; then
  AC_MSG_ERROR([the specified C++ compiler does not appear to support namespaces])
fi

AX_CXX_HAVE_STL
if test "$ax_cv_cxx_have_stl" != "yes" ; then
  AC_MSG_ERROR([the specified C++ compiler does not appear support the Standard Template Library])
fi

AC_MSG_CHECKING([whether compiler supports __builtin_expect keyword])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
if (__builtin_expect(true,1)) return 0;
]])],[
have_builtin_expect=1
AC_MSG_RESULT(yes)],[
have_builtin_expect=0
AC_MSG_RESULT(no)])
AC_DEFINE_UNQUOTED(HAVE_BUILTIN_EXPECT,$have_builtin_expect,[Boolean value indicating whether the C++ compiler supports the __builtin_expect keyword])

# Checks for other programs.
AC_PROG_RANLIB

AC_ARG_WITH([ar],
  AS_HELP_STRING([--with-ar=ARG],[manually set ar to ARG]))
AC_CHECK_PROGS(AR, $with_ar ar,,)
if test -z "$AR" ; then
  AC_MSG_ERROR([ar not found, specify via --with-ar=ARG])
fi

AC_ARG_WITH([M4],
  AS_HELP_STRING([--with-M4=ARG],[manually set M4 to ARG]))
AC_CHECK_PROGS(M4, $with_M4 gm4 m4,,)
if test -z "$M4" ; then
  AC_MSG_ERROR([m4 preprocessor not found, specify via --with-M4=ARG])
fi

AC_ARG_WITH([doxygen-dir],
  AS_HELP_STRING([--with-doxygen-dir=PATH],[manually set directory where the doxygen executable resides to PATH]),
  DOXYGEN_PATH=$withval
  AC_CHECK_PROGS(DOXYGEN, "$DOXYGEN_PATH/doxygen" doxygen,,),
  DOXYGEN_PATH=""
  AC_CHECK_PROGS(DOXYGEN, doxygen,,))
if test -z "$DOXYGEN" ; then
  AC_MSG_WARN([doxygen not found])
  echo "if doxygen is installed, specify its location via --with-doxygen-dir=PATH"
  HAVE_DOXYGEN="NO"
else
  HAVE_DOXYGEN="YES"
fi
AC_SUBST(HAVE_DOXYGEN,$HAVE_DOXYGEN)
AC_SUBST(DOXYGEN,$DOXYGEN)
AC_SUBST(DOXYGEN_PATH,$DOXYGEN_PATH)

AC_ARG_WITH([dot-dir],
  AS_HELP_STRING([--with-dot-dir=PATH],[manually set directory where the Graphviz dot executable resides to PATH]),
  DOT_PATH=$withval
  AC_CHECK_PROGS(DOT, "$DOT_PATH/dot" dot,,),
  DOT_PATH=""
  AC_CHECK_PROGS(DOT, dot,,))
if test -z "$DOT" ; then
  AC_MSG_WARN([dot not found])
  echo "if the Graphviz dot tool is installed, specify its location via --with-dot-dir=PATH"
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

# Checks to see whether subdomain indices should be enabled.
AC_ARG_ENABLE([subdomain-indices],
  AS_HELP_STRING(--enable-subdomain-indices,[enable the use of subdomain indices for standard IB mechanical elements @<:@default=no@:>@]),
  [SUBDOMAIN_INDICES=$enable_subdomain_indices], [SUBDOMAIN_INDICES=no])

if test "$SUBDOMAIN_INDICES" == "yes" ; then
  echo "subdomain indices are ENABLED!"
  AC_DEFINE([ENABLE_SUBDOMAIN_INDICES],[true],[Boolean value indicating whether we are using subdomain indices])
else
  echo "subdomain indices are DISABLED!"
  AC_DEFINE([ENABLE_SUBDOMAIN_INDICES],[false],[Boolean value indicating whether we are using subdomain indices])
fi

###########################################################################
# Checks for libraries.
###########################################################################

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

AC_ARG_WITH([hdf5],
  AS_HELP_STRING(--with-hdf5=PATH,location of required hdf5 installation @<:@default=/usr/local/sfw/hdf5/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-hdf5=PATH])
   fi
   HDF5_DIR=$withval],
  [HDF5_DIR=/usr/local/sfw/hdf5/current])

AC_ARG_WITH([hypre],
  AS_HELP_STRING(--with-hypre=PATH,location of required hypre installation @<:@default=/usr/local/sfw/hypre/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-hypre=PATH])
   fi
   HYPRE_DIR=$withval],
  [HYPRE_DIR=/usr/local/sfw/hypre/current])

AC_ARG_ENABLE([libmesh],
  AS_HELP_STRING(--enable-libmesh,enable support for the optional libMesh library @<:@default=no@:>@),
  [LIBMESH_ENABLED=$enableval], [LIBMESH_ENABLED=no])
if test "$LIBMESH_ENABLED" == "" ; then
  LIBMESH_ENABLED=no
fi

AC_ARG_WITH([libmesh],
  AS_HELP_STRING([--with-libmesh=PATH],[location of optional libMesh installation @<:@default=/usr/local/sfw/libmesh/current@:>@]),
  [if test "$LIBMESH_ENABLED" == "no" ; then
     AC_MSG_WARN([--with-libmesh is specified, but support for libMesh is DISABLED])
     LIBMESH_DIR=NONE
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([you must specify an existing directory when using --with-libmesh=PATH])
     fi
     LIBMESH_DIR=$withval
   fi],
  [LIBMESH_DIR=/usr/local/sfw/libmesh/current])

AC_ARG_WITH([libmesh-compile-mode],
  AS_HELP_STRING([--with-libmesh-compile-mode=METHOD],[libMesh compile mode options include opt, devel, dbg, prof, oprof @<:@default=dbg@:>@]),
  [if test "$LIBMESH_ENABLED" == "no" ; then
     AC_MSG_WARN([--with-libmesh-compile-mode is specified, but support for libMesh is DISABLED])
     LIBMESH_COMPILE_MODE=NONE
   else
     LIBMESH_COMPILE_MODE=$withval
   fi],
  [LIBMESH_COMPILE_MODE=dbg])

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

AC_ARG_WITH([samrai],
  AS_HELP_STRING(--with-samrai=PATH,location of required SAMRAI installation @<:@default=/usr/local/sfw/samrai/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-samrai=PATH])
   fi
   SAMRAI_DIR=$withval],
  [SAMRAI_DIR=/usr/local/sfw/samrai/current])

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

###########################################################################
# Configure libraries.
###########################################################################

echo
echo "========================================"
echo "Configuring required package BLAS/LAPACK"
echo "========================================"
if (test -r $PETSC_DIR/$PETSC_ARCH/conf/petscvariables) ; then
  PETSC_BLASLAPACK_LIB=`grep "BLASLAPACK_LIB" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | sed -e 's/.*BLASLAPACK_LIB = //'`
  if test "x$PETSC_BLASLAPACK_LIB" != x ; then
    if test "x$with_blas" != "x" -a "$with_blas" != "$PETSC_BLASLAPACK_LIB" ; then
      AC_MSG_WARN([BLAS flags specified by --with-blas appear to disagree with PETSc BLAS/LAPACK flags; using PETSc settings])
    fi
    if test "x$with_lapack" != "x" -a "$with_lapack" != "$PETSC_BLASLAPACK_LIB" ; then
      AC_MSG_WARN([LAPACK flags specified by --with-lapack appear to disagree with PETSc BLAS/LAPACK flags; using PETSc settings])
    fi
    echo "using PETSc BLAS/LAPACK library flags: $PETSC_BLASLAPACK_LIB"
    with_blas="$PETSC_BLASLAPACK_LIB"
    with_lapack=""
  fi
  PETSC_BLASLAPACK_INCLUDE=`grep "BLASLAPACK_INCLUDE" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | sed -e 's/.*BLASLAPACK_INCLUDE = //'`
  if test "x$PETSC_BLASLAPACK_INCLUDE" != x ; then
    CPPFLAGS="${PETSC_BLASLAPACK_INCLUDE} $CPPFLAGS"
  fi
fi
AX_BLAS(,)
AX_LAPACK(,)
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FCLIBS"

AM_CONDITIONAL([BLITZ_ENABLED],[test "$BLITZ_ENABLED" == "yes"])
if test "$BLITZ_ENABLED" == "yes"; then
  echo
  echo "===================================="
  echo "Configuring optional package Blitz++"
  echo "===================================="
  CONFIGURE_BLITZ
else
  echo
  echo "===================================="
  echo "Optional package Blitz++ is DISABLED"
  echo "===================================="
fi

AM_CONDITIONAL([SILO_ENABLED],[test "$SILO_ENABLED" == "yes"])
if test "$SILO_ENABLED" == "yes"; then
  echo
  echo "================================="
  echo "Configuring optional package Silo"
  echo "================================="
  CONFIGURE_SILO
else
  echo
  echo "================================="
  echo "Optional package Silo is DISABLED"
  echo "================================="
fi

echo
echo "================================="
echo "Configuring required package hdf5"
echo "================================="
CONFIGURE_HDF5

echo
echo "=================================="
echo "Configuring required package hypre"
echo "=================================="
CONFIGURE_HYPRE

echo
echo "=================================="
echo "Configuring required package PETSc"
echo "=================================="
CONFIGURE_PETSC

AM_CONDITIONAL([LIBMESH_ENABLED],[test "$LIBMESH_ENABLED" == "yes"])
if test "$LIBMESH_ENABLED" == "yes"; then
  echo
  echo "===================================="
  echo "Configuring optional package libMesh"
  echo "===================================="
  CONFIGURE_LIBMESH
else
  echo
  echo "===================================="
  echo "Optional package libMesh is DISABLED"
  echo "===================================="
fi

echo
echo "==================================="
echo "Configuring required package SAMRAI"
echo "==================================="
CONFIGURE_SAMRAI

###########################################################################
# Output files.
###########################################################################
echo
echo "================"
echo "Outputting files"
echo "================"
AC_CONFIG_FILES([
  Makefile
  bin/Makefile
  config/make.inc
  doc/Doxyfile
  doc/Makefile
  doc/online_Doxyfile
  examples/Makefile
  examples/IB/Makefile
  examples/IB/cell_centered/Makefile
  examples/IB/cell_centered/ex0/Makefile
  examples/IB/cell_centered/ex1/Makefile
  examples/IB/cell_centered/ex2/Makefile
  examples/IB/staggered_grid/Makefile
  examples/IB/staggered_grid/explicit/Makefile
  examples/IB/staggered_grid/explicit/ex0/Makefile
  examples/IB/staggered_grid/explicit/ex1/Makefile
  examples/IB/staggered_grid/explicit/ex2/Makefile
  examples/IB/staggered_grid/explicit/ex3/Makefile
  examples/IB/staggered_grid/explicit/ex4/Makefile
  examples/IB/staggered_grid/implicit/Makefile
  examples/IB/staggered_grid/implicit/ex0/Makefile
  examples/IB/staggered_grid/implicit/ex1/Makefile
  examples/IBFE/Makefile
  examples/IBFE/explicit/Makefile
  examples/IBFE/explicit/ex0/Makefile
  examples/IBFE/explicit/ex1/Makefile
  examples/IBFE/explicit/ex2/Makefile
  examples/IBFE/explicit/ex3/Makefile
  examples/IBFE/explicit/ex4/Makefile
  examples/IBFE/explicit/ex5/Makefile
  examples/IBFE/explicit/ex6/Makefile
  examples/adv_diff/Makefile
  examples/advect/Makefile
  examples/navier_stokes/Makefile
  examples/navier_stokes/cell_centered/Makefile
  examples/navier_stokes/cell_centered/approx_projection/Makefile
  examples/navier_stokes/cell_centered/approx_projection/ex0/Makefile
  examples/navier_stokes/cell_centered/approx_projection/ex1/Makefile
  examples/navier_stokes/cell_centered/approx_projection/ex2/Makefile
  examples/navier_stokes/cell_centered/exact_projection/Makefile
  examples/navier_stokes/cell_centered/exact_projection/ex0/Makefile
  examples/navier_stokes/cell_centered/exact_projection/ex1/Makefile
  examples/navier_stokes/staggered_grid/Makefile
  examples/navier_stokes/staggered_grid/ex0/Makefile
  examples/navier_stokes/staggered_grid/ex1/Makefile
  examples/navier_stokes/staggered_grid/ex2/Makefile
  examples/navier_stokes/staggered_grid/ex3/Makefile
  examples/vc_navier_stokes/Makefile
  examples/vc_navier_stokes/operator_test/Makefile
  examples/vc_navier_stokes/solver_test/Makefile
  lib/Makefile
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
  src/tools/Makefile
])
AC_CONFIG_SUBDIRS([ibtk])
AC_OUTPUT
echo
echo "===================================================================="
echo "===================================================================="
echo
AC_MSG_NOTICE([IBAMR appears to have been configured successfully.  Now:

   make lib        Build the IBAMR library and its supporting IBTK library
   make check      Build and run basic library tests
   make examples   Build all of the IBAMR and IBTK example programs

Note that the IBAMR and IBTK example programs are not run by "make examples".
To run the examples, e.g.,

   make examples
   cd examples/IB/cell_centered/ex0
   ./main2d input2d

You may also build and run individual example programs, e.g., via

   cd examples/IB/cell_centered/ex0
   make examples
   ./main2d input2d])
