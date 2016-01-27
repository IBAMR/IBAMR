###########################################################################
# Prelude.
###########################################################################
AC_INIT([IBAMR],[140227],[griffith@cims.nyu.edu],[IBAMR],[https://ibamr.googlecode.com])
AC_PREREQ([2.69])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config/IBAMR_config.h.tmp])
AX_PREFIX_CONFIG_H([config/IBAMR_config.h],[IBAMR],[config/IBAMR_config.h.tmp])
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

###########################################################################
# Checks for optional and required third-party libraries.
###########################################################################
PACKAGE_INITIALIZE_ENVIRONMENT
CONFIGURE_LIBMESH # must configure libMesh first to scope its configuration
CONFIGURE_BOOST("$ABSOLUTE_SRCDIR/ibtk/contrib/boost","\$(abs_top_builddir)/ibtk/contrib/boost")
CONFIGURE_EIGEN("$ABSOLUTE_SRCDIR/ibtk/contrib/eigen","\$(abs_top_builddir)/ibtk/contrib/eigen")
CONFIGURE_HDF5
CONFIGURE_HYPRE
CONFIGURE_MUPARSER("$ABSOLUTE_SRCDIR/ibtk/contrib/muparser","\$(abs_top_builddir)/ibtk/contrib/muparser")
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
  config/make.inc
  examples/Makefile
  examples/ConstraintIB/Makefile
  examples/ConstraintIB/eel2d/Makefile
  examples/ConstraintIB/eel3d/Makefile  
  examples/IB/Makefile
  examples/IB/explicit/Makefile
  examples/IB/explicit/ex0/Makefile
  examples/IB/explicit/ex1/Makefile
  examples/IB/explicit/ex2/Makefile
  examples/IB/explicit/ex3/Makefile
  examples/IB/explicit/ex4/Makefile
  examples/IB/explicit/ex5/Makefile
  examples/IB/explicit/ex6/Makefile
  examples/IB/implicit/Makefile
  examples/IB/implicit/ex0/Makefile
  examples/IB/implicit/ex1/Makefile
  examples/IB/implicit/ex2/Makefile
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
  examples/IMP/Makefile
  examples/IMP/explicit/Makefile
  examples/IMP/explicit/ex0/Makefile
  examples/adv_diff/Makefile
  examples/adv_diff/ex0/Makefile
  examples/adv_diff/ex1/Makefile
  examples/adv_diff/ex2/Makefile
  examples/advect/Makefile
  examples/navier_stokes/Makefile
  examples/navier_stokes/ex0/Makefile
  examples/navier_stokes/ex1/Makefile
  examples/navier_stokes/ex2/Makefile
  examples/navier_stokes/ex3/Makefile
  examples/navier_stokes/ex4/Makefile
  examples/navier_stokes/ex5/Makefile
  examples/navier_stokes/ex6/Makefile
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
  tests/Makefile
  tests/Stokes/Makefile
  tests/Stokes/test0/Makefile
  tests/Stokes-IB/Makefile
  tests/Stokes-IB/test0/Makefile
  tests/Stokes-IB/test1/Makefile
  tests/Stokes-IB/test2/Makefile
])
AC_CONFIG_SUBDIRS([ibtk])
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
