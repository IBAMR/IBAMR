# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_HYPRE],[
echo
echo "=================================="
echo "Configuring required package HYPRE"
echo "=================================="

PACKAGE_SETUP_ENVIRONMENT

if test `grep -c HYPRE "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables"` != 0 ; then
  AC_MSG_NOTICE([PETSc appears to provide hypre; using PETSc hypre library])
  PETSC_BUNDLES_HYPRE=yes
else
  PETSC_BUNDLES_HYPRE=no

  AC_ARG_WITH([hypre],
    AS_HELP_STRING(--with-hypre=PATH,location of required hypre installation),
    [if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-hypre=PATH])
     fi
     HYPRE_DIR=$withval])

  if test x$HYPRE_DIR != x ; then
    if test -d "${HYPRE_DIR}/include" ; then
      HYPRE_CPPFLAGS="-I${HYPRE_DIR}/include"
    fi
    if test -d "${HYPRE_DIR}/lib" ; then
      HYPRE_LDFLAGS="-L${HYPRE_DIR}/lib"
    fi
  fi

  CPPFLAGS_PREPEND($HYPRE_CPPFLAGS)
  AC_CHECK_HEADER([HYPRE.h],,AC_MSG_ERROR([could not find header file HYPRE.h]))

  LDFLAGS_PREPEND($HYPRE_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([HYPRE])
  if test "$HAVE_LIBHYPRE" = no ; then
    AC_MSG_ERROR([could not find working libHYPRE])
  fi

  PACKAGE_CPPFLAGS_PREPEND($HYPRE_CPPFLAGS)
  PACKAGE_LDFLAGS_PREPEND($HYPRE_LDFLAGS)
  PACKAGE_LIBS_PREPEND("$LIBHYPRE")
fi

PACKAGE_RESTORE_ENVIRONMENT

])
