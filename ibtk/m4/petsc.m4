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
#include <petscversion.h>
]], [[
#if (PETSC_VERSION_(3,3,0))
#else
asdf
#endif
]])],[PETSC_VERSION_3_3_0=yes],[PETSC_VERSION_3_3_0=no])

if test "$PETSC_VERSION_3_3_0" == "yes"; then
  echo "using PETSc 3.3.0"
else
  AC_MSG_ERROR([incorrect PETSc version detected: please use PETSc 3.3.0])
fi

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <petsc.h>
]], [[
#ifdef PETSC_HAVE_X
#if (PETSC_HAVE_X == 1)
asdf
#endif
#endif
]])],[PETSC_HAVE_X=no],[PETSC_HAVE_X=yes])

if test "$PETSC_HAVE_X" == "yes"; then
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
])