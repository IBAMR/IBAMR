## ---------------------------------------------------------------------
##
## Copyright (c) 2011 - 2022 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

###########################################################################
# Prelude.
###########################################################################
AC_INIT([IBTK],[171104],[boyceg@gmail.com],[IBAMR],[https://ibamr.github.io])
AC_PREREQ([2.69])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config/IBTK_config.h.tmp])
AX_PREFIX_CONFIG_H([config/IBTK_config.h],[IBTK],[config/IBTK_config.h.tmp])
# Also copy our own CMake compatibility shim header:
AC_CONFIG_LINKS([config/ibtk/config.h:$srcdir/config/config.h.tmp])
AM_INIT_AUTOMAKE([1.12 -Wall -Werror dist-bzip2 foreign -Wno-extra-portability subdir-objects])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([disable])
AC_LANG(C++)
dnl libMesh will automatically detect and use C++14 features if they are
dnl available by default (e.g., with GCC 6.0 and newer). Ensure that we have the
dnl same version flag by using exactly the same macro setup as libMesh:
dnl
dnl (this variable is set up by libMesh)
enablecxx11=yes
AX_CXX_COMPILE_STDCXX_11([],[mandatory])

###########################################################################
# Checks for programs.
###########################################################################
AC_ARG_VAR(LDFLAGS,[linker flags, e.g. -L<lib dir> if you have libraries in a nonstandard directory <lib dir>])
AC_ARG_VAR(LIBS,[linker commands, e.g. -l<lib> if you require a nonstandard library <lib>])
AC_ARG_VAR(FCLIBS,[linker flags and commands, e.g. -L<lib dir> and -l<lib> for the "Fortran intrinsic and run-time libraries" that are required to successfully link a Fortran program or shared library])
AC_PROG_CC		# Check for a C compiler.
AC_PROG_CPP		# Check for a C preprocessor.
AC_PROG_CXX		# Check for a C++ compiler.
AC_PROG_CXXCPP		# Check for a C++ preprocessor.
AC_PROG_FC		# Check for a Fortran compiler.
AC_FC_LIBRARY_LDFLAGS	# Determine the linker flags for the Fortran intrinsic and runtime libraries.
AC_FC_WRAPPERS		# Determine the form of the symbol name mangling used by the Fortran compiler and setup wrappers to perform the name mangling.
F77="$FC"		# Set F77, FFLAGS, and FLIBS to avoid needing to call AC_PROG_F77.
FFLAGS="$FCFLAGS"
FLIBS="$FCLIBS"
AC_FC_SRCEXT(f)
AC_SUBST(F77)
AC_SUBST(FFLAGS)
AC_SUBST(FLIBS)
AX_PROG_CC_MPI		# Check for an MPI C compiler wrapper and use it as the default C compiler.
AX_PROG_CXX_MPI		# Check for an MPI C++ compiler wrapper and use it as the default C++ compiler.
AC_HEADER_MAJOR
AC_HEADER_STDC
LT_INIT
AC_LIB_PROG_LD_GNU
AC_ARG_WITH([M4],
  AS_HELP_STRING([--with-M4=ARG],[manually set M4 to ARG]))
AC_CHECK_PROGS(M4, $with_M4 gm4 m4,,)
if test -z "$M4" ; then
  AC_MSG_ERROR([m4 preprocessor not found, specify via --with-M4=ARG])
fi
AC_PROG_GREP
AC_PROG_SED
CHECK_BUILTIN_EXPECT
CHECK_BUILTIN_PREFETCH
CHECK_PRAGMA_KEYWORD
CONFIGURE_DOXYGEN
CONFIGURE_DOT

###########################################################################
# Checks for optional and required third-party libraries.
###########################################################################
PACKAGE_INITIALIZE_ENVIRONMENT
echo
echo "=================================="
echo "Configuring required package Boost"
echo "=================================="
AC_MSG_NOTICE([first check for a system Boost library; if not found, revert to bundled Boost library])
BOOST_REQUIRE([1.66.0],[
  AC_MSG_WARN([could not find system Boost library, using bundled Boost library])
  USING_BUNDLED_BOOST=yes])

# Configure libraries that manage other dependencies of IBAMR:
dnl libMesh must be first since it will override PETSc settings (in the case
dnl that libMesh uses one version of PETSc and another version of PETSc was
dnl provided by argument)
CONFIGURE_LIBMESH
CONFIGURE_PETSC
# Configure other libraries that are not themselves dependencies:
CONFIGURE_BOOST("$ABSOLUTE_SRCDIR/contrib/boost","\$(abs_top_builddir)/contrib/boost")
CONFIGURE_EIGEN("$ABSOLUTE_SRCDIR/contrib/eigen","\$(abs_top_builddir)/contrib/eigen")
CONFIGURE_HDF5
CONFIGURE_HYPRE
CONFIGURE_MUPARSER("$ABSOLUTE_SRCDIR/contrib/muparser","\$(abs_top_builddir)/contrib/muparser")
# configure dependencies of dependencies:
CONFIGURE_SILO
CONFIGURE_SAMRAI
PACKAGE_SETUP_ENVIRONMENT
LIBS="$LIBS $PACKAGE_CONTRIB_LIBS"

# Check that SAMRAI and LIBMESH have mutually compatible debug settings:
if test "$LIBMESH_ENABLED" = yes; then
    if test "$SAMRAI_DEBUG_CHECK_ASSERTIONS" = true; then
        if test "$LIBMESH_METHOD" = dbg || test "$LIBMESH_METHOD" = devel; then
        : # do nothing
        else
            AC_MSG_ERROR([The detected SAMRAI library was compiled with debugging support enabled but the detected libMesh library was not. This incompatibility is not allowed since different debug macros can lead to linkage errors. Please reconfigure IBAMR with either a devel or debug libMesh method.])
        fi
        else # SAMRAI is in release mode
        if test "$LIBMESH_METHOD" = dbg || test "$LIBMESH_METHOD" = devel; then
            AC_MSG_ERROR([The detected SAMRAI library was compiled without debugging support enabled but the detected libMesh library was. This incompatibility is not allowed since different debug macros can lead to linkage errors. Please reconfigure IBAMR with either an opt, pro, or oprof libMesh method.])
        fi
    fi
fi


###########################################################################
# Additional header configuration.
###########################################################################

dnl We have to use a funny name for the key to guarantee that this is at the top.
AH_VERBATIM(A_INCLUDE_MACRO_HEADER, [/* Warn that this file is deprecated. */
#ifndef IBTK_SKIP_CONFIG_DEPRECATION_WARNINGS
#pragma message("This file has been deprecated: use ibtk/config.h instead")
#endif
])

dnl We have to use a funny name for the key to guarantee that it comes after
dnl all constant definitions.
AH_VERBATIM(Z_INCLUDE_MACRO_HEADER, [/* Include the new configuration file. */
#include <ibtk/config.h>])

###########################################################################
# Output files.
###########################################################################
echo
echo "================"
echo "Outputting files"
echo "================"

AC_CONFIG_FILES([
Makefile
examples/Makefile
examples/CCLaplace/Makefile
examples/CCPoisson/Makefile
examples/GhostCells/Makefile
examples/PETScOps/Makefile
examples/PETScOps/ProlongationMat/Makefile
examples/PhysBdryOps/Makefile
examples/SCLaplace/Makefile
examples/SCPoisson/Makefile
examples/VCLaplace/Makefile
examples/VCViscousSolver/Makefile
lib/Makefile
scripts/Makefile
src/Makefile
src/boundary/Makefile
src/boundary/cf_interface/Makefile
src/boundary/cf_interface/fortran/Makefile
src/boundary/physical_boundary/Makefile
src/boundary/physical_boundary/fortran/Makefile
src/coarsen_ops/Makefile
src/coarsen_ops/fortran/Makefile
src/fortran/Makefile
src/lagrangian/Makefile
src/lagrangian/fortran/Makefile
src/math/Makefile
src/math/fortran/Makefile
src/refine_ops/Makefile
src/refine_ops/fortran/Makefile
src/solvers/Makefile
src/solvers/impls/Makefile
src/solvers/impls/fortran/Makefile
src/solvers/interfaces/Makefile
src/solvers/wrappers/Makefile
src/utilities/Makefile
src/utilities/fortran/Makefile
contrib/Makefile
contrib/muparser/Makefile
contrib/muparser/src/Makefile
])

AC_OUTPUT
