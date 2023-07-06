## ---------------------------------------------------------------------
##
## Copyright (c) 2006 - 2022 by the IBAMR developers
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
AC_INIT([IBAMR],[171104],[boyceg@gmail.com],[IBAMR],[https://ibamr.github.io])
AC_PREREQ([2.69])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config/IBAMR_config.h.tmp])
AX_PREFIX_CONFIG_H([config/IBAMR_config.h],[IBAMR],[config/IBAMR_config.h.tmp])
# Also copy our own CMake compatibility shim header:
AC_CONFIG_LINKS([config/ibamr/config.h:$srcdir/config/config.h.tmp])
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

###########################################################################
# Version information (requires sed).
###########################################################################
dnl Since this is an M4 script we must use the quadrigraphs @<:@ and @:>@ for
dnl brackets, i.e., for [ and ]
ibamr_version_major=$($SED -E "s/^(@<:@0-9@:>@+)\..*/\1/" "$srcdir/VERSION")
ibamr_version_minor=$($SED -E "s/^@<:@^.@:>@+\.(@<:@0-9@:>@+)\.@<:@0-9@:>@+.*/\1/" "$srcdir/VERSION")
ibamr_version_subminor=$($SED -E "s/^@<:@^.@:>@+\.@<:@^.@:>@+\.(@<:@0-9@:>@+).*/\1/" "$srcdir/VERSION")
AC_DEFINE_UNQUOTED(VERSION_MAJOR, $ibamr_version_major, [Major version of IBAMR.])
AC_DEFINE_UNQUOTED(VERSION_MINOR, $ibamr_version_minor, [Minor version of IBAMR.])
AC_DEFINE_UNQUOTED(VERSION_SUBMINOR, $ibamr_version_subminor, [Subminor version of IBAMR.])

dnl We have to use a funny name for the key to guarantee that this is at the top.
AH_VERBATIM(A_INCLUDE_MACRO_HEADER, [/* Warn that this file is deprecated. */
#ifndef IBAMR_SKIP_CONFIG_DEPRECATION_WARNINGS
#pragma message("This file has been deprecated: use ibamr/config.h instead")
#endif
])

dnl We have to use a funny name for the key to guarantee that it comes after
dnl all constant definitions.
AH_VERBATIM(Z_INCLUDE_MACRO_HEADER, [/* Include the new configuration file. */
#include <ibamr/config.h>])

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
CONFIGURE_BOOST("$ABSOLUTE_SRCDIR/ibtk/contrib/boost","\$(abs_top_builddir)/ibtk/contrib/boost")
CONFIGURE_EIGEN("$ABSOLUTE_SRCDIR/ibtk/contrib/eigen","\$(abs_top_builddir)/ibtk/contrib/eigen")
CONFIGURE_HDF5
CONFIGURE_HYPRE
CONFIGURE_MUPARSER("$ABSOLUTE_SRCDIR/ibtk/contrib/muparser","\$(abs_top_builddir)/ibtk/contrib/muparser")
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

# old version:
AC_SUBST(MPIEXEC)
AC_SUBST(NUMDIFF)

# for compatibility with CMake, define new versions
MPIEXEC_EXECUTABLE="$MPIEXEC"
NUMDIFF_EXECUTABLE="$NUMDIFF"
AC_SUBST(MPIEXEC_EXECUTABLE)
AC_SUBST(NUMDIFF_EXECUTABLE)

###########################################################################
# Output files.
###########################################################################
echo
echo "================"
echo "Outputting files"
echo "================"
AC_CONFIG_FILES([
  Makefile
  config/make.inc
  examples/Makefile
  examples/CIB/Makefile
  examples/CIB/ex0/Makefile
  examples/CIB/ex1/Makefile
  examples/CIB/ex2/Makefile
  examples/CIB/ex3/Makefile
  examples/CIB/ex4/Makefile
  examples/ConstraintIB/Makefile
  examples/ConstraintIB/eel2d/Makefile
  examples/ConstraintIB/eel3d/Makefile
  examples/ConstraintIB/falling_sphere/Makefile
  examples/ConstraintIB/flow_past_cylinder/Makefile
  examples/ConstraintIB/impulsively_started_cylinder/Makefile
  examples/ConstraintIB/knifefish/Makefile
  examples/ConstraintIB/moving_plate/Makefile
  examples/ConstraintIB/oscillating_rigid_cylinder/Makefile
  examples/ConstraintIB/stokes_first_problem/Makefile
  examples/DLM/Makefile
  examples/DLM/ex0/Makefile
  examples/IB/Makefile
  examples/IB/explicit/Makefile
  examples/IB/explicit/ex0/Makefile
  examples/IB/explicit/ex1/Makefile
  examples/IB/explicit/ex2/Makefile
  examples/IB/explicit/ex3/Makefile
  examples/IB/explicit/ex4/Makefile
  examples/IB/explicit/ex5/Makefile
  examples/IB/explicit/ex6/Makefile
  examples/IB/explicit/ex7/Makefile
  examples/IBFE/Makefile
  examples/IBFE/explicit/Makefile
  examples/IBFE/explicit/ex0/Makefile
  examples/IBFE/explicit/ex1/Makefile
  examples/IBFE/explicit/ex2/Makefile
  examples/IBFE/explicit/ex3/Makefile
  examples/IBFE/explicit/ex4/Makefile
  examples/IBFE/explicit/ex5/Makefile
  examples/IBFE/explicit/ex6/Makefile
  examples/IBFE/explicit/ex7/Makefile
  examples/IBFE/explicit/ex8/Makefile
  examples/IBFE/explicit/ex9/Makefile
  examples/IBFE/explicit/ex10/Makefile
  examples/IBFE/explicit/ex11/Makefile
  examples/IBLevelSet/Makefile
  examples/IBLevelSet/ex0/Makefile
  examples/IBLevelSet/ex1/Makefile
  examples/IBLevelSet/ex2/Makefile
  examples/IBLevelSet/ex3/Makefile
  examples/IBLevelSet/ex4/Makefile
  examples/IBLevelSet/ex5/Makefile
  examples/IIM/Makefile
  examples/IIM/ex0/Makefile
  examples/IIM/ex1/Makefile
  examples/IIM/ex2/Makefile
  examples/IIM/ex3/Makefile
  examples/IIM/ex4/Makefile
  examples/IIM/ex5/Makefile
  examples/IIM/ex6/Makefile
  examples/IIM/ex7/Makefile
  examples/IIM/ex8/Makefile
  examples/IIM/ex9/Makefile
  examples/IMP/Makefile
  examples/IMP/explicit/Makefile
  examples/IMP/explicit/ex0/Makefile
  examples/adv_diff/Makefile
  examples/adv_diff/ex0/Makefile
  examples/adv_diff/ex1/Makefile
  examples/adv_diff/ex2/Makefile
  examples/adv_diff/ex3/Makefile
  examples/adv_diff/ex4/Makefile
  examples/adv_diff/ex5/Makefile
  examples/adv_diff/ex6/Makefile
  examples/adv_diff/ex7/Makefile
  examples/advect/Makefile
  examples/complex_fluids/Makefile
  examples/complex_fluids/ex0/Makefile
  examples/complex_fluids/ex1/Makefile
  examples/complex_fluids/ex2/Makefile
  examples/complex_fluids/ex3/Makefile
  examples/complex_fluids/ex4/Makefile
  examples/fe_mechanics/Makefile
  examples/fe_mechanics/ex0/Makefile
  examples/level_set/Makefile
  examples/level_set/ex0/Makefile
  examples/level_set/ex1/Makefile
  examples/multiphase_flow/Makefile
  examples/multiphase_flow/ex0/Makefile
  examples/multiphase_flow/ex1/Makefile
  examples/multiphase_flow/ex2/Makefile
  examples/multiphase_flow/ex3/Makefile
  examples/multiphase_flow/ex4/Makefile
  examples/multiphase_flow/ex5/Makefile
  examples/multiphase_flow/ex6/Makefile
  examples/multiphase_flow/ex7/Makefile
  examples/multiphase_flow/ex8/Makefile
  examples/multiphase_flow/ex9/Makefile
  examples/multiphase_flow/ex10/Makefile
  examples/multiphase_flow/ex11/Makefile
  examples/multiphase_flow/ex12/Makefile
  examples/multiphase_flow/ex13/Makefile
  examples/navier_stokes/Makefile
  examples/navier_stokes/ex0/Makefile
  examples/navier_stokes/ex1/Makefile
  examples/navier_stokes/ex2/Makefile
  examples/navier_stokes/ex3/Makefile
  examples/navier_stokes/ex4/Makefile
  examples/navier_stokes/ex5/Makefile
  examples/navier_stokes/ex6/Makefile
  examples/vc_navier_stokes/Makefile
  examples/vc_navier_stokes/ex0/Makefile
  examples/vc_navier_stokes/ex1/Makefile
  examples/vc_navier_stokes/ex2/Makefile
  examples/wave_tank/Makefile
  examples/wave_tank/ex0/Makefile
  examples/wave_tank/ex1/Makefile
  lib/Makefile
  src/Makefile
  src/fortran/Makefile
  src/IB/Makefile
  src/adv_diff/Makefile
  src/adv_diff/fortran/Makefile
  src/advect/Makefile
  src/advect/fortran/Makefile
  src/complex_fluids/Makefile
  src/complex_fluids/fortran/Makefile
  src/level_set/Makefile
  src/level_set/fortran/Makefile
  src/navier_stokes/Makefile
  src/navier_stokes/fortran/Makefile
  src/phase_change/Makefile
  src/utilities/Makefile
  src/wave_generation/Makefile
  tests/attest.conf
  tests/Makefile
  tests/adv_diff/Makefile
  tests/advect/Makefile
  tests/coarsen/Makefile
  tests/complex_fluids/Makefile
  tests/CIB/Makefile
  tests/ConstraintIB/Makefile
  tests/fe_mechanics/Makefile
  tests/IB/Makefile
  tests/IBFE/Makefile
  tests/IIM/Makefile
  tests/IBTK/Makefile
  tests/IMP/Makefile
  tests/interpolate/Makefile
  tests/level_set/Makefile
  tests/multiphase_flow/Makefile
  tests/navier_stokes/Makefile
  tests/physical_boundary/Makefile
  tests/refine/Makefile
  tests/spread/Makefile
  tests/vc_navier_stokes/Makefile
  tests/wave_tank/Makefile
])
AC_CONFIG_SUBDIRS([ibtk])

ibamr_version_string=$(cat "$srcdir/VERSION")
AC_MSG_NOTICE([Configured version of IBAMR is $ibamr_version_string])
# check if git is available:
which git>/dev/null
FOUND_GIT=$?
if test "$FOUND_GIT" -eq 0 # 0 indicates 'which git' succeeded
then
# and we have a git repository:
    if test -f $srcdir/.git/HEAD
    then
        AC_MSG_NOTICE([IBAMR configured with git commit $(git rev-parse HEAD)])
    else
        AC_MSG_NOTICE([Current source directory is not a git repository. Skipping git version info.])
    fi
else
    AC_MSG_NOTICE([Current source directory is not a git repository. Skipping git version info.])
fi

AC_OUTPUT
echo
echo "===================================================================="
echo "===================================================================="
echo
AC_MSG_NOTICE([IBAMR appears to have been configured successfully.  Now:

   make lib        Build the IBAMR library and its supporting IBTK library
   make examples   Build all of the IBAMR and IBTK example programs

Note that the IBAMR and IBTK example programs are not run by "make examples".
To run the examples, e.g.,

   make examples
   cd examples/IB/explicit/ex0
   ./main2d input2d

You may also build and run individual example programs, e.g., via

   cd examples/IB/explicit/ex0
   make examples
   ./main2d input2d])
