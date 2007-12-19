###########################################################################
# Prelude.
###########################################################################
AC_INIT([IBAMR],[070517],[griffith AT cims DOT nyu DOT edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config/IBAMR_config.h])
AM_INIT_AUTOMAKE([-Wall -Werror dist-bzip2 foreign])

AC_LANG(C++)		# Set the default language to C++.

###########################################################################
# Checks for programs.
###########################################################################
AC_ARG_VAR(LIBS,[linker commands, e.g. -l<lib> if you require a nonstandard library <lib>])
AC_ARG_VAR(FLIBS,[linker flags and commands, e.g. -L<lib dir> and -l<lib> for the "Fortran 77 intrinsic and run-time libraries" that are required to successfully link a Fortran 77 program or shared library])
AC_PROG_CC		# Check for a C compiler.
AC_PROG_CXX		# Check for a C++ compiler.
AC_PROG_F77		# Check for a Fortran compiler.
AC_F77_LIBRARY_LDFLAGS	# Determine the linker flags for the Fortran intrinsic and runtime libraries.
AC_F77_WRAPPERS

# Check rpath syntax.
AC_ARG_WITH([rpath],
  AS_HELP_STRING([--with-rpath],[indicate whether additional linker flags, e.g. -rpath<lib dir> are required when linking with shared libraries]),
  [USE_RPATH_FLAG=$withval],
  [USE_RPATH_FLAG="yes"])

AC_ARG_WITH([rpath-syntax],
  AS_HELP_STRING([--with-rpath-syntax=ARG],[linker flags, e.g. -rpath<lib dir> are required when linking with shared libraries @<:@default=-Wl,-rpath,@:>@]),
  [RPATH_FLAG=$withval],
  [RPATH_FLAG="-Wl,-rpath,"])

if test "$USE_RPATH_FLAG" == "yes" ; then
  echo "using rpath flags when linking against shared libraries"
  echo "rpath flag syntax is: ${RPATH_FLAG}<lib dir>"
else
  echo "not using rpath flags when linking against shared libraries"
fi

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
ACX_MPI(,AC_MSG_WARN([unable to find working MPI installation]))
CC="$MPICC"
AC_LANG_POP([C])

ACX_MPI(,AC_MSG_WARN([unable to find working MPI installation]))
CXX="$MPICXX"

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
###AC_DISABLE_SHARED
###AC_PROG_LIBTOOL

AC_ARG_WITH([M4],
  AS_HELP_STRING([--with-M4=ARG],[manually set M4 to ARG]))
AC_CHECK_PROGS(M4, $with_M4 gm4 m4,,)
if test -z "$M4" ; then
  AC_MSG_ERROR([m4 preprocessor not found, specify via --with-M4=ARG])
fi

AC_ARG_WITH([doxygen-dir],
  AS_HELP_STRING([--with-doxygen-dir],[manually set directory where the doxygen executable resides to ARG]),
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
  echo "  if you wish to enable support for large ghost cell widths, specify --enable-large-ghost-cell-width=yes"
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
  echo "  if you wish to enable support for large ghost cell widths, specify --enable-expensive-cf-interpolation=yes"
  AC_DEFINE([USING_EXPENSIVE_CF_INTERPOLATION],[false],[Boolean value indicating whether we are using the expensive version of the quadratic coarse-fine interface code])
fi

###########################################################################
# Checks for libraries.
###########################################################################

# Required external libraries.

# Check for BLAS and LAPACK libraries.
ACX_BLAS(,)
ACX_LAPACK(,)
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

# Check for zlib library and configure zlib options.
echo
echo "================================="
echo "Configuring required package zlib"
echo "================================="
AC_ARG_WITH([zlib],
  AS_HELP_STRING(--with-zlib=PATH,location of required zlib installation @<:@default=/usr/local/sfw/zlib/current@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-zlib=PATH])
   fi
   ZLIB_DIR=$withval],
  [ZLIB_DIR=/usr/local/sfw/zlib/current])

if test -d "${ZLIB_DIR}/lib" ; then
  if test "${USE_RPATH_FLAG}" == "yes" ; then
    LDFLAGS="${RPATH_FLAG}${ZLIB_DIR}/lib $LDFLAGS"
  fi
  LDFLAGS="-L${ZLIB_DIR}/lib $LDFLAGS"
fi
if test -d "${ZLIB_DIR}/include" ; then
  CPPFLAGS="-I${ZLIB_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([zlib.h],,AC_MSG_WARN([could not find header file zlib.h]))

AC_LANG_PUSH([C])
AC_CHECK_LIB([z],[inflateEnd],,AC_MSG_WARN([could not find libz]))
AC_LANG_POP([C])

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
  if test "${USE_RPATH_FLAG}" == "yes" ; then
    LDFLAGS="${RPATH_FLAG}${HDF5_DIR}/lib $LDFLAGS"
  fi
  LDFLAGS="-L${HDF5_DIR}/lib $LDFLAGS"
fi
if test -d "${HDF5_DIR}/include" ; then
  CPPFLAGS="-I${HDF5_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([hdf5.h],,AC_MSG_WARN([could not find header file hdf5.h]))

AC_LANG_PUSH([C])
AC_CHECK_LIB([hdf5],[H5open],,
  echo "checking to see if we need to include szip to use hdf5..."
  AC_CHECK_LIB([sz], [SZ_encoder_enabled])
  unset ac_cv_lib_hdf5_H5open
  AC_CHECK_LIB([hdf5],[H5open],,AC_MSG_WARN([could not find libhdf5])))
AC_CHECK_LIB([hdf5_hl],[H5open],,AC_MSG_WARN([could not find libhdf5_hl]))
AC_LANG_POP([C])

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
  if test "${USE_RPATH_FLAG}" == "yes" ; then
    LDFLAGS="${RPATH_FLAG}${HYPRE_DIR}/lib $LDFLAGS"
  fi
  LDFLAGS="-L${HYPRE_DIR}/lib $LDFLAGS"
fi
if test -d "${HYPRE_DIR}/include" ; then
  CPPFLAGS="-I${HYPRE_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([HYPRE.h],,AC_MSG_WARN([could not find header file HYPRE.h]))

AC_LANG_PUSH([C])
AC_CHECK_LIB([HYPRE_utilities],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_utilities]))
AC_CHECK_LIB([HYPRE_multivector],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_multivector]))
AC_CHECK_LIB([HYPRE_krylov],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_krylov]))
AC_CHECK_LIB([HYPRE_seq_mv],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_seq_mv]))
AC_CHECK_LIB([HYPRE_parcsr_mv],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_parcsr_mv]))
AC_CHECK_LIB([HYPRE_IJ_mv],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_IJ_mv]))
AC_CHECK_LIB([HYPRE_DistributedMatrix],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_DistributedMatrix]))
AC_CHECK_LIB([HYPRE_MatrixMatrix],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_MatrixMatrix]))
AC_CHECK_LIB([HYPRE_Euclid],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_Euclid]))
AC_CHECK_LIB([HYPRE_ParaSails],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_ParaSails]))
AC_CHECK_LIB([HYPRE_DistributedMatrixPilutSolver],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_DistributedMatrixPilutSolver]))
AC_CHECK_LIB([HYPRE_parcsr_block_mv],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_parcsr_block_mv]))
AC_CHECK_LIB([HYPRE_parcsr_ls],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_parcsr_ls]))
AC_CHECK_LIB([HYPRE_struct_mv],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_struct_mv]))
AC_CHECK_LIB([HYPRE_struct_ls],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_struct_ls]))
AC_CHECK_LIB([HYPRE_sstruct_mv],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_sstruct_mv]))
AC_CHECK_LIB([HYPRE_sstruct_ls],[hypre_CAlloc],,AC_MSG_WARN([could not find libHYPRE_sstruct_ls]))
AC_LANG_POP([C])

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
  if test "${USE_RPATH_FLAG}" == "yes" ; then
    LDFLAGS="${RPATH_FLAG}${SAMRAI_DIR}/lib $LDFLAGS"
  fi
  LDFLAGS="-L${SAMRAI_DIR}/lib $LDFLAGS"
fi
if test -d "${SAMRAI_DIR}/include" ; then
  CPPFLAGS="-I${SAMRAI_DIR}/include $CPPFLAGS"
fi

AC_SUBST(SAMRAI_DIR,[${SAMRAI_DIR}])
AC_SUBST(SAMRAI_FORTDIR,[${SAMRAI_DIR}/include])

AC_CHECK_HEADER([SAMRAI_config.h],,AC_MSG_WARN([could not find header file SAMRAI_config.h]))
AC_CHECK_LIB([SAMRAI],[main])

AC_ARG_ENABLE([samrai-2d],
  AS_HELP_STRING(--enable-samrai-2d,enable optional support for two-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI2D=$enablevar], [SAMRAI2D=yes])
AM_CONDITIONAL([SAMRAI2D_ENABLED],[test "$SAMRAI2D" == "yes"])

if test "$SAMRAI2D" == "yes"; then
  AC_CHECK_LIB([SAMRAI2d_hier],[main],,AC_MSG_WARN([could not find libSAMRAI2d_hier]))
  AC_CHECK_LIB([SAMRAI2d_xfer],[main],,AC_MSG_WARN([could not find libSAMRAI2d_xfer]))
  AC_CHECK_LIB([SAMRAI2d_pdat_special],[main],,AC_MSG_WARN([could not find libSAMRAI2d_pdat_special]))
  AC_CHECK_LIB([SAMRAI2d_pdat_std],[main],,AC_MSG_WARN([could not find libSAMRAI2d_pdat_std]))
  AC_CHECK_LIB([SAMRAI2d_math_special],[main],,AC_MSG_WARN([could not find libSAMRAI2d_math_special]))
  AC_CHECK_LIB([SAMRAI2d_math_std],[main],,AC_MSG_WARN([could not find libSAMRAI2d_math_std]))
  AC_CHECK_LIB([SAMRAI2d_mesh],[main],,AC_MSG_WARN([could not find libSAMRAI2d_mesh]))
  AC_CHECK_LIB([SAMRAI2d_mblk],[main],,AC_MSG_WARN([could not find libSAMRAI2d_mblk]))
  AC_CHECK_LIB([SAMRAI2d_geom],[main],,AC_MSG_WARN([could not find libSAMRAI2d_geom]))
  AC_CHECK_LIB([SAMRAI2d_solv],[main],,AC_MSG_WARN([could not find libSAMRAI2d_solv]))
  AC_CHECK_LIB([SAMRAI2d_algs],[main],,AC_MSG_WARN([could not find libSAMRAI2d_algs]))
  AC_CHECK_LIB([SAMRAI2d_appu],[main],,AC_MSG_WARN([could not find libSAMRAI2d_appu]))
else
  AC_MSG_WARN([configuring without the two-dimensional SAMRAI library])
fi

AC_ARG_ENABLE([samrai-3d],
  AS_HELP_STRING(--enable-samrai-3d,enable optional support for three-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI3D=$enablevar], [SAMRAI3D=yes])
AM_CONDITIONAL([SAMRAI3D_ENABLED],[test "$SAMRAI3D" == "yes"])

if test "$SAMRAI3D" == "yes"; then
  AC_CHECK_LIB([SAMRAI3d_hier],[main],,AC_MSG_WARN([could not find libSAMRAI3d_hier]))
  AC_CHECK_LIB([SAMRAI3d_xfer],[main],,AC_MSG_WARN([could not find libSAMRAI3d_xfer]))
  AC_CHECK_LIB([SAMRAI3d_pdat_special],[main],,AC_MSG_WARN([could not find libSAMRAI3d_pdat_special]))
  AC_CHECK_LIB([SAMRAI3d_pdat_std],[main],,AC_MSG_WARN([could not find libSAMRAI3d_pdat_std]))
  AC_CHECK_LIB([SAMRAI3d_math_special],[main],,AC_MSG_WARN([could not find libSAMRAI3d_math_special]))
  AC_CHECK_LIB([SAMRAI3d_math_std],[main],,AC_MSG_WARN([could not find libSAMRAI3d_math_std]))
  AC_CHECK_LIB([SAMRAI3d_mesh],[main],,AC_MSG_WARN([could not find libSAMRAI3d_mesh]))
  AC_CHECK_LIB([SAMRAI3d_mblk],[main],,AC_MSG_WARN([could not find libSAMRAI3d_mblk]))
  AC_CHECK_LIB([SAMRAI3d_geom],[main],,AC_MSG_WARN([could not find libSAMRAI3d_geom]))
  AC_CHECK_LIB([SAMRAI3d_solv],[main],,AC_MSG_WARN([could not find libSAMRAI3d_solv]))
  AC_CHECK_LIB([SAMRAI3d_algs],[main],,AC_MSG_WARN([could not find libSAMRAI3d_algs]))
  AC_CHECK_LIB([SAMRAI3d_appu],[main],,AC_MSG_WARN([could not find libSAMRAI3d_appu]))
else
  AC_MSG_WARN([configuring without the three-dimensional SAMRAI library])
fi

# Determine how to use the RobinBcCoefStrategy class.

AC_COMPILE_IFELSE(
  [
#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <RobinBcCoefStrategy.h>
#include <Patch.h>
#include <Variable.h>

#define NDIM 2

using namespace SAMRAI;
using namespace SAMRAI::hier;
using namespace SAMRAI::pdat;
using namespace SAMRAI::solv;
using namespace SAMRAI::tbox;

void foo()
{
    Pointer<ArrayData<NDIM,double> > acoef_data = NULL;
    Pointer<ArrayData<NDIM,double> > bcoef_data = NULL;
    Pointer<ArrayData<NDIM,double> > gcoef_data = NULL;
    Pointer<Variable<NDIM> > var = NULL;
    Box<NDIM> patch_box;
    Patch<NDIM> patch(patch_box, NULL);
    BoundaryBox<NDIM> bdry_box;
    double fill_time = 0.0;

    Pointer<RobinBcCoefStrategy<NDIM> > bc_coef;
    bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, var,
        patch, bdry_box, fill_time);

    return;
}
  ],
  [AC_DEFINE([USING_OLD_ROBIN_BC_INTERFACE],[0],
             [Flag indicating whether we are using the RobinBcCoefStrategy interface from SAMRAI 2.1])],
  [AC_DEFINE([USING_OLD_ROBIN_BC_INTERFACE],[1],
             [Flag indicating whether we are using the RobinBcCoefStrategy interface from SAMRAI 2.1])])

# Check for PETSc library and configure PETSc options.
echo
echo "=================================="
echo "Configuring required package PETSc"
echo "=================================="
AC_ARG_VAR(PETSC_DIR,[the location of the PETSc installation that is to be used.  Note that multiple versions of PETSc may be installed within the same filesystem, with each version corresponding to a different value of PETSC_DIR.])
AC_ARG_VAR(PETSC_ARCH,[the PETSc configuration that is to be used, corresponding to the configuration files located in ${PETSC_DIR}/bmake/${PETSC_ARCH}.  Note that multiple configurations of PETSc can be installed for a particular version of the PETSc library.  Each configuration will correspond to a different value of PETSC_ARCH.])
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
echo "using PETSc lib     directory = ${PETSC_DIR}/lib/${PETSC_ARCH}"
echo "using PETSc include directory = ${PETSC_DIR}/include, ${PETSC_DIR}/bmake/${PETSC_ARCH}"

if test -d "${PETSC_DIR}/lib/${PETSC_ARCH}" ; then
  if test "${USE_RPATH_FLAG}" == "yes" ; then
    LDFLAGS="${RPATH_FLAG}${PETSC_DIR}/lib/${PETSC_ARCH} $LDFLAGS"
  fi
  LDFLAGS="-L${PETSC_DIR}/lib/${PETSC_ARCH} $LDFLAGS"
fi
if test -d "${PETSC_DIR}/include" ; then
  CPPFLAGS="-I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/${PETSC_ARCH} $CPPFLAGS"
fi

AC_CHECK_HEADER([petsc.h],,AC_MSG_WARN([could not find header file petsc.h]))

AC_LANG_PUSH([C])
AC_CHECK_LIB([petsc],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetsc])
              missing_petsc_lib=yes])
AC_CHECK_LIB([petscvec],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetscvec])
              missing_petsc_lib=yes])
AC_CHECK_LIB([petscmat],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetscmat])
              missing_petsc_lib=yes])
AC_CHECK_LIB([petscdm],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetscdm])
              missing_petsc_lib=yes])
AC_CHECK_LIB([petscksp],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetscksp])
              missing_petsc_lib=yes])
AC_CHECK_LIB([petscsnes],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetscsnes])
              missing_petsc_lib=yes])
AC_CHECK_LIB([petscts],[PetscTrMalloc],,
             [AC_MSG_WARN([could not find libpetscts])
              missing_petsc_lib=yes])
AC_LANG_POP([C])

if test "$missing_petsc_lib" == "yes"; then
  echo "could not find one or more of the PETSc libraries"
  echo "trying to link using the C++ compiler instead of the C compiler"

  unset ac_cv_lib_petsc_PetscTrMalloc
  unset ac_cv_lib_petscvec_PetscTrMalloc
  unset ac_cv_lib_petscmat_PetscTrMalloc
  unset ac_cv_lib_petscdm_PetscTrMalloc
  unset ac_cv_lib_petscksp_PetscTrMalloc
  unset ac_cv_lib_petscsnes_PetscTrMalloc
  unset ac_cv_lib_petscts_PetscTrMalloc

  AC_CHECK_LIB([petsc],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetsc])
                missing_petsc_lib=yes])
  AC_CHECK_LIB([petscvec],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetscvec])
                missing_petsc_lib=yes])
  AC_CHECK_LIB([petscmat],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetscmat])
                missing_petsc_lib=yes])
  AC_CHECK_LIB([petscdm],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetscdm])
                missing_petsc_lib=yes])
  AC_CHECK_LIB([petscksp],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetscksp])
                missing_petsc_lib=yes])
  AC_CHECK_LIB([petscsnes],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetscsnes])
                missing_petsc_lib=yes])
  AC_CHECK_LIB([petscts],[PetscTrMalloc],,
               [AC_MSG_WARN([could not find libpetscts])
                missing_petsc_lib=yes])
fi

# Optional libraries.

# Check for Blitz++ library and configure Blitz++ options.
AC_ARG_ENABLE([blitz],
  AS_HELP_STRING(--enable-blitz,enable support for the optional Blitz++ library @<:@default=yes@:>@),
  [BLITZ_ENABLED=$enableval], [BLITZ_ENABLED=yes])
if test "$BLITZ_ENABLED" == "" ; then
  BLITZ_ENABLED=yes
fi

AC_ARG_WITH([blitz],
  AS_HELP_STRING(--with-blitz=PATH,location of optional Blitz++ installation @<:@default=/usr/local/sfw/blitz/current@:>@),
  [if test "$BLITZ_ENABLED" = no ; then
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
    if test "${USE_RPATH_FLAG}" == "yes" ; then
      LDFLAGS="${RPATH_FLAG}${BLITZ_DIR}/lib $LDFLAGS"
    fi
    LDFLAGS="-L${BLITZ_DIR}/lib $LDFLAGS"
  fi
  if test -d "${BLITZ_DIR}/include" ; then
   CPPFLAGS="-I${BLITZ_DIR}/include $CPPFLAGS"
  fi

  AC_CHECK_HEADER([blitz/blitz.h],,AC_MSG_WARN([could not find header file blitz/blitz.h]))
  AC_CHECK_LIB([blitz],[main],,AC_MSG_WARN([could not find libblitz]))
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
  [if test "$SILO_ENABLED" = no ; then
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
    ###if test "${USE_RPATH_FLAG}" == "yes" ; then
    ###  LDFLAGS="${RPATH_FLAG}${SILO_DIR}/lib $LDFLAGS"
    ###fi
    LDFLAGS="-L${SILO_DIR}/lib $LDFLAGS"
  fi
  if test -d "${SILO_DIR}/include" ; then
    CPPFLAGS="-I${SILO_DIR}/include $CPPFLAGS"
  fi

  AC_CHECK_HEADER([silo.h],,AC_MSG_WARN([could not find header file silo.h]))

  AC_LANG_PUSH([C])
  AC_CHECK_LIB([silo],[DBVersion],,AC_MSG_WARN([could not find libsilo]))
  AC_LANG_POP([C])
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
AC_CONFIG_SUBDIRS([stools])
AC_CONFIG_FILES([
  Makefile
  config/make.inc
  doc/Makefile doc/Doxyfile
  sample/IB/Makefile sample/empty_project/Makefile
  src/Makefile src/fortran/Makefile
  src/IB/Makefile
  src/adv_diff/Makefile src/adv_diff/fortran/Makefile
  src/advect/Makefile src/advect/fortran/Makefile
  src/lagrangian/Makefile src/lagrangian/fortran/Makefile
  src/navier_stokes/Makefile src/navier_stokes/fortran/Makefile
  src/utilities/Makefile src/utilities/get_version.C
  src/tests/Makefile
  src/tests/IB/Makefile
  src/tests/IB_multilevel/Makefile
  src/tests/IB_target/Makefile
  src/tests/adv_diff/Makefile
  src/tests/advect/Makefile
  src/tests/navier_stokes/Makefile
  src/tests/navier_stokes/analytic_soln/Makefile
  src/tests/navier_stokes/lid_driven_cavity/Makefile
  src/tests/navier_stokes/periodic_shear_layer/Makefile
  src/tests/pIB/Makefile
  src/tools/Makefile
])
AC_OUTPUT
