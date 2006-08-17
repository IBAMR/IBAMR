###########################################################################
# Prelude.
###########################################################################
AC_INIT([IBAMR],[0.1],[griffith AT cims DOT nyu DOT edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_HEADERS([config/IBAMR_config.h])
AM_INIT_AUTOMAKE([-Wall -Werror dist-bzip2 foreign])

AC_LANG(C++)		# Set the default language to C++.

###########################################################################
# Checks for programs.
###########################################################################
AC_PROG_CC		# Check for a C compiler.
AC_PROG_CXX		# Check for a C++ compiler.
AC_PROG_F77		# Check for a Fortran compiler.
AC_F77_LIBRARY_LDFLAGS	# Determine the linker flags for the Fortran intrinsic and runtime libraries.
AC_F77_WRAPPERS

# Check for MPICXX and configure MPI options.
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
AC_ARG_WITH([mpicxx],
  AS_HELP_STRING(--with-mpicxx=ARG,MPI C++ compiler command),
  [MPICXX=$withval],)

ACX_MPI(,AC_MSG_WARN([unable to find working MPI installation]))
CXX="$MPICXX"
if test "$MPILIBS" != "" ; then
  LIBS="$MPILIBS $LIBS"
fi

# Checks for other programs.
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

AC_ARG_WITH([M4],
   AS_HELP_STRING([--with-M4=ARG],[manually set M4 to ARG]))
AC_CHECK_PROGS(M4, $with_M4 gm4 m4,,)
if test -z "$M4" ; then
  AC_MSG_ERROR([m4 preprocessor not found, specify via --with-M4=ARG])
fi

###########################################################################
# Checks for libraries.
###########################################################################

# Required external libraries.

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
  LDFLAGS="-L${SAMRAI_DIR}/lib -R${SAMRAI_DIR}/lib $LDFLAGS"
fi
if test -d "${SAMRAI_DIR}/include" ; then
  CPPFLAGS="-I${SAMRAI_DIR}/include $CPPFLAGS"
fi

AC_SUBST(SAMRAI_FORTDIR,[${SAMRAI_DIR}/include])

AC_CHECK_HEADER([SAMRAI_config.h],,AC_MSG_WARN([could not find header file SAMRAI_config.h]))
AC_CHECK_LIB([SAMRAI],[malloc])

AC_ARG_ENABLE([samrai-1d],
  AS_HELP_STRING(--enable-samrai-1d,enable optional support for one-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI1D=$enableval], [SAMRAI1D=yes])

if test "$SAMRAI1D" = yes; then
  AC_CHECK_LIB([SAMRAI1d_hier],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_hier]))
  AC_CHECK_LIB([SAMRAI1d_xfer],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_xfer]))
  AC_CHECK_LIB([SAMRAI1d_pdat_special],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_pdat_special]))
  AC_CHECK_LIB([SAMRAI1d_pdat_std],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_pdat_std]))
  AC_CHECK_LIB([SAMRAI1d_math_special],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_math_special]))
  AC_CHECK_LIB([SAMRAI1d_math_std],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_math_std]))
  AC_CHECK_LIB([SAMRAI1d_mesh],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_mesh]))
  AC_CHECK_LIB([SAMRAI1d_mblk],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_mblk]))
  AC_CHECK_LIB([SAMRAI1d_geom],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_geom]))
  AC_CHECK_LIB([SAMRAI1d_solv],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_solv]))
  AC_CHECK_LIB([SAMRAI1d_algs],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_algs]))
  AC_CHECK_LIB([SAMRAI1d_appu],[malloc],,AC_MSG_WARN([could not find libSAMRAI1d_appu]))
else
  AC_MSG_WARN([configuring without the one-dimensional SAMRAI library])
fi

AC_ARG_ENABLE([samrai-2d],
  AS_HELP_STRING(--enable-samrai-2d,enable optional support for two-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI2D=$enablevar], [SAMRAI2D=yes])

if test "$SAMRAI2D" = yes; then
  AC_CHECK_LIB([SAMRAI2d_hier],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_hier]))
  AC_CHECK_LIB([SAMRAI2d_xfer],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_xfer]))
  AC_CHECK_LIB([SAMRAI2d_pdat_special],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_pdat_special]))
  AC_CHECK_LIB([SAMRAI2d_pdat_std],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_pdat_std]))
  AC_CHECK_LIB([SAMRAI2d_math_special],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_math_special]))
  AC_CHECK_LIB([SAMRAI2d_math_std],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_math_std]))
  AC_CHECK_LIB([SAMRAI2d_mesh],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_mesh]))
  AC_CHECK_LIB([SAMRAI2d_mblk],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_mblk]))
  AC_CHECK_LIB([SAMRAI2d_geom],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_geom]))
  AC_CHECK_LIB([SAMRAI2d_solv],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_solv]))
  AC_CHECK_LIB([SAMRAI2d_algs],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_algs]))
  AC_CHECK_LIB([SAMRAI2d_appu],[malloc],,AC_MSG_WARN([could not find libSAMRAI2d_appu]))
else
  AC_MSG_WARN([configuring without the two-dimensional SAMRAI library])
fi

AC_ARG_ENABLE([samrai-3d],
  AS_HELP_STRING(--enable-samrai-3D,enable optional support for three-dimensional SAMRAI objects @<:@default=yes@:>@),
  [SAMRAI3D=$enablevar], [SAMRAI3D=yes])

if test "$SAMRAI3D" = yes; then
  AC_CHECK_LIB([SAMRAI3d_hier],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_hier]))
  AC_CHECK_LIB([SAMRAI3d_xfer],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_xfer]))
  AC_CHECK_LIB([SAMRAI3d_pdat_special],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_pdat_special]))
  AC_CHECK_LIB([SAMRAI3d_pdat_std],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_pdat_std]))
  AC_CHECK_LIB([SAMRAI3d_math_special],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_math_special]))
  AC_CHECK_LIB([SAMRAI3d_math_std],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_math_std]))
  AC_CHECK_LIB([SAMRAI3d_mesh],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_mesh]))
  AC_CHECK_LIB([SAMRAI3d_mblk],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_mblk]))
  AC_CHECK_LIB([SAMRAI3d_geom],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_geom]))
  AC_CHECK_LIB([SAMRAI3d_solv],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_solv]))
  AC_CHECK_LIB([SAMRAI3d_algs],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_algs]))
  AC_CHECK_LIB([SAMRAI3d_appu],[malloc],,AC_MSG_WARN([could not find libSAMRAI3d_appu]))
else
  AC_MSG_WARN([configuring without the three-dimensional SAMRAI library])
fi

# Check for PETSc library and configure PETSc options.
echo
echo "=================================="
echo "Configuring required package PETSc"
echo "=================================="
AC_ARG_VAR(PETSC_DIR,[location of the PETSc installation to be used])
AC_ARG_VAR(PETSC_ARCH,[name of the PETSc configuration to be used])
AC_ARG_WITH([petsc],
  AS_HELP_STRING(--with-petsc=PATH,location of required PETSc installation @<:@default=PETSC_DIR@:>@),
  [if test ! -d "$withval" ; then
     AC_MSG_ERROR([you must specify an existing directory when using --with-petsc=PATH])
   fi
   PETSC_DIR=$withval],)
AC_ARG_WITH([PETSC_ARCH],
  AS_HELP_STRING(--with-petsc-arch=ARCH,PETSc architecture @<:@default=PETSC_ARCH@:>@),
  [PETSC_ARCH=$withval],)
AC_SUBST(PETSC_DIR,$PETSC_DIR)
AC_SUBST(PETSC_ARCH,$PETSC_ARCH)

echo "using PETSc architecture = ${PETSC_ARCH}"
echo "using PETSc root    directory = ${PETSC_DIR}"
echo "using PETSc lib     directory = ${PETSC_DIR}/lib/${PETSC_ARCH}"
echo "using PETSc include directory = ${PETSC_DIR}/include, ${PETSC_DIR}/bmake/${PETSC_ARCH}"

if test -d "${PETSC_DIR}/lib/${PETSC_ARCH}" ; then
  LDFLAGS="-L${PETSC_DIR}/lib/${PETSC_ARCH} -R${PETSC_DIR}/lib/${PETSC_ARCH} $LDFLAGS"
fi
if test -d "${PETSC_DIR}/include" ; then
  CPPFLAGS="-I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/${PETSC_ARCH} $CPPFLAGS"
fi

AC_CHECK_HEADER([petsc.h],,AC_MSG_WARN([could not find header file petsc.h]))
AC_CHECK_LIB([petsc],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetsc]))
AC_CHECK_LIB([petscvec],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetscvec]))
AC_CHECK_LIB([petscmat],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetscmat]))
AC_CHECK_LIB([petscdm],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetscdm]))
AC_CHECK_LIB([petscksp],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetscksp]))
AC_CHECK_LIB([petscsnes],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetscsnes]))
AC_CHECK_LIB([petscts],[PetscTrMalloc],,AC_MSG_WARN([could not find libpetscts]))

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
  LDFLAGS="-L${HYPRE_DIR}/lib -R${HYPRE_DIR}/lib $LDFLAGS"
fi
if test -d "${HYPRE_DIR}/include" ; then
  CPPFLAGS="-I${HYPRE_DIR}/include $CPPFLAGS"
fi

AC_CHECK_HEADER([HYPRE.h],,AC_MSG_WARN([could not find header file HYPRE.h]))
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
  LDFLAGS="-L${HDF5_DIR}/lib -R${HDF5_DIR}/lib $LDFLAGS"
fi
if test -d "${HDF5_DIR}/include" ; then
  CPPFLAGS="-I${HDF5_DIR}/include $CPPFLAGS"  
fi

AC_CHECK_HEADER([hdf5.h],,AC_MSG_WARN([could not find header file hdf5.h]))
AC_CHECK_LIB([hdf5],[H5open],,AC_MSG_WARN([could not find libhdf5]))

# Optional libraries.

# Check for Blitz++ library and configure Blitz++ options.
AC_ARG_ENABLE([blitz],
  AS_HELP_STRING(--enable-blitz,enable support for the optional Blitz++ library @<:@default=yes@:>@),
  [BLITZ_ENABLED=$enableval], [BLITZ_ENABLED=yes])
if test "$BLITZ_ENABLED" = "" ; then
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

if test "$BLITZ_ENABLED" = yes; then
  echo
  echo "===================================="
  echo "Configuring optional package Blitz++"
  echo "===================================="

  if test -d "${BLITZ_DIR}/lib" ; then
    LDFLAGS="-L${BLITZ_DIR}/lib -R${BLITZ_DIR}/lib $LDFLAGS"
  fi
  if test -d "${BLITZ_DIR}/include" ; then
   CPPFLAGS="-I${BLITZ_DIR}/include $CPPFLAGS"
  fi

  AC_CHECK_HEADER([blitz/blitz.h],,AC_MSG_WARN([could not find header file blitz/blitz.h]))
  AC_CHECK_LIB([blitz],[_blitz_id],,AC_MSG_WARN([could not find libblitz]))
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
if test "$SILO_ENABLED" = "" ; then
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

if test "$SILO_ENABLED" = yes; then
  echo
  echo "================================="
  echo "Configuring optional package Silo"
  echo "================================="

  if test -d "${SILO_DIR}/lib" ; then
    LDFLAGS="-L${SILO_DIR}/lib -R${SILO_DIR}/lib $LDFLAGS"
  fi
  if test -d "${SILO_DIR}/include" ; then
    CPPFLAGS="-I${SILO_DIR}/include $CPPFLAGS"  
  fi

  AC_CHECK_HEADER([silo.h],,AC_MSG_WARN([could not find header file silo.h]))
  AC_CHECK_LIB([silo],[DBVersion],,AC_MSG_WARN([could not find libsilo]))
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
AC_CONFIG_FILES([
  Makefile
  src/Makefile src/fortran/Makefile
  src/advect/Makefile src/advect/fortran/Makefile src/advect/tests/Makefile
  src/utilities/Makefile
])
AC_OUTPUT
