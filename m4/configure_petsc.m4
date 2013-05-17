# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PETSC],[
echo
echo "=================================="
echo "Configuring required package PETSc"
echo "=================================="

PACKAGE_SETUP_ENVIRONMENT

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

AC_MSG_NOTICE([using PETSC_DIR: ${PETSC_ARCH}])
AC_MSG_NOTICE([using PETSc_ARCH: ${PETSC_ARCH}])

PETSC_CC_INCLUDES=`grep "PETSC_CC_INCLUDES =" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | sed -e 's/.*=//' | sed -e 's/^[ \t]*//'`
PETSC_EXTERNAL_LIB_BASIC=`grep "PETSC_EXTERNAL_LIB_BASIC =" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | sed -e 's/.*=//' | sed -e 's/^[ \t]*//'`
PETSC_WITH_EXTERNAL_LIB=`grep "PETSC_WITH_EXTERNAL_LIB =" $PETSC_DIR/$PETSC_ARCH/conf/petscvariables | sed -e 's/.*=//' | sed -e 's/^[ \t]*//'`

CPPFLAGS_PREPEND($PETSC_CC_INCLUDES)
AC_CHECK_HEADER([petsc.h],,AC_MSG_ERROR([could not find header file petsc.h]))

AC_MSG_CHECKING([for PETSc version 3.3])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <petscversion.h>
]], [[
#if (PETSC_VERSION_(3,3,0))
#else
asdf
#endif
]])],[PETSC_VERSION_3_3_0=yes],[PETSC_VERSION_3_3_0=no])
AC_MSG_RESULT([PETSC_VERSION_3_3_0])
if test "$PETSC_VERSION_3_3_0" = no; then
  AC_MSG_ERROR([incorrect PETSc version detected: please use PETSc 3.3])
fi

LIBS_PREPEND($PETSC_EXTERNAL_LIB_BASIC)
AC_LIB_HAVE_LINKFLAGS([petsc])
if test "$HAVE_LIBPETSC" = no ; then
  AC_MSG_ERROR([could not find working libpetsc])
fi

PACKAGE_CPPFLAGS_PREPEND($PETSC_CC_INCLUDES)
PACKAGE_LIBS_PREPEND($PETSC_WITH_EXTERNAL_LIB)

PACKAGE_RESTORE_ENVIRONMENT

])