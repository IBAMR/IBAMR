###########################################################################
# Prelude.
###########################################################################
AC_INIT([IBTK],[140227],[griffith@cims.nyu.edu],[IBTK],[https://ibamr.googlecode.com])
AC_PREREQ([2.69])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config/IBTK_config.h.tmp])
AX_PREFIX_CONFIG_H([config/IBTK_config.h],[IBTK],[config/IBTK_config.h.tmp])
AM_INIT_AUTOMAKE([1.12 -Wall -Werror dist-bzip2 foreign -Wno-extra-portability subdir-objects])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([disable])
AC_LANG(C++)

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
CONFIGURE_DOXYGEN
CONFIGURE_DOT

###########################################################################
# Checks for optional and required third-party libraries.
###########################################################################
PACKAGE_INITIALIZE_ENVIRONMENT
CONFIGURE_LIBMESH # must configure libMesh first to scope its configuration
CONFIGURE_BOOST("$ABSOLUTE_SRCDIR/contrib/boost","\$(abs_top_builddir)/contrib/boost")
CONFIGURE_EIGEN("$ABSOLUTE_SRCDIR/contrib/eigen","\$(abs_top_builddir)/contrib/eigen")
CONFIGURE_HDF5
CONFIGURE_HYPRE
CONFIGURE_MUPARSER("$ABSOLUTE_SRCDIR/contrib/muparser","\$(abs_top_builddir)/contrib/muparser")
CONFIGURE_PETSC
CONFIGURE_SAMRAI
CONFIGURE_SILO
CONFIGURE_GSL
PACKAGE_SETUP_ENVIRONMENT
LIBS="$LIBS $PACKAGE_CONTRIB_LIBS"

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
  examples/PETScOps/Makefile
  examples/PETScOps/ProlongationMat/Makefile
  examples/PhysBdryOps/Makefile
  examples/SCLaplace/Makefile
  examples/SCPoisson/Makefile
  examples/VCLaplace/Makefile
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
  contrib/Makefile
])
if test "$USING_BUNDLED_MUPARSER" = yes ; then
  AC_CONFIG_SUBDIRS([contrib/muparser])
fi
AC_OUTPUT
