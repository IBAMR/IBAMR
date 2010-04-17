AC_DEFUN([CONFIGURE_PETSC],[

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

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <petsc.h>
]], [[
#ifdef PETSC_HAVE_X11
#if (PETSC_HAVE_X11 == 1)
asdf
#endif
#endif
]])],[PETSC_HAS_X11=no],[PETSC_HAS_X11=yes])

if test "$PETSC_HAS_X11" == "yes"; then
  echo "NOTE: PETSc appears to have been configured with X11 support enabled"
  echo "looking for libX11"
  LDFLAGS="-L/usr/X11/lib $LDFLAGS"
  AC_LIB_HAVE_LINKFLAGS([X11])
  if test "$HAVE_LIBX11" != "yes"; then
    AC_MSG_ERROR([could not find working X11 library
try modifying the value of the LDFLAGS environment variable,
or reconfigure PETSc with the flag --with-x=0])
  fi
  LIBS="$LIBX11 $LIBS"
fi

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
])